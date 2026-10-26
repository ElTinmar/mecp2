import os
import re
from multiprocessing import Pool
import cv2
import numpy as np
import pandas as pd
from tqdm import tqdm


# Parameters
FRAME_BUFFER_INTERVAL = 5  # Frames
WINDOW_SEC = 15 * 60 + 20
MIN_PREY_SIZE = 20  # Pixels
MAX_PREY_SIZE = 750  # Pixels

# Paths
ROOT_DIR = '.'
INPUT_FOLDER = 'raw_videos'
OUTPUT_FOLDER = 'processed_videos'

input_path = os.path.join(ROOT_DIR, INPUT_FOLDER)
output_path = os.path.join(ROOT_DIR, OUTPUT_FOLDER)
os.makedirs(output_path, exist_ok=True)


# --- Load existing processed files if CSV exists ---
existing_filenames = set()
video_data_csv = os.path.join(output_path, 'video_data.csv')
counts_csv = os.path.join(output_path, 'prey_counts.csv')

if os.path.exists(video_data_csv):
    existing_video_df = pd.read_csv(video_data_csv)
    existing_filenames.update(existing_video_df['file_name'].tolist())
else:
    existing_video_df = pd.DataFrame()

if os.path.exists(counts_csv):
    existing_count_df = pd.read_csv(counts_csv)
else:
    existing_count_df = pd.DataFrame()


# Functions
def process_video_wrapper(args):
    return process_video(*args)

def process_video(video_path, output_dir, delay=None):
    """
    Process a single video for a fixed 15-minute window after a delay.

    Args:
        video_path (str): Path to the input .avi video.
        output_dir (str): Directory to store the output detection video.
        delay (int | float): Delayed start time in SECONDS.

    Returns:
        filename (str)
        all_prey_counts (list of int)
    """

    # Extract filename and prepare output video path
    filename = os.path.basename(video_path)
    output_file = os.path.join(
        output_dir, filename.replace(".avi", "_detection.avi")
    )

    # Open video and get FPS
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        print(f"Failed to open video: {video_path}")
        return filename, []

    fps = cap.get(cv2.CAP_PROP_FPS)
    if fps <= 0:
        print(f"Invalid FPS for video: {video_path}")
        cap.release()
        return filename, []

    # -----------------------------
    # Frame window definition
    # -----------------------------
    delay = delay or 0  # seconds

    start_frame = int(delay * fps)
    window_frames = int(WINDOW_SEC * fps)
    end_frame = start_frame + window_frames

    total_frames = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))

    # SAFEGUARD: video too short
    if end_frame > total_frames:
        print(
            f"Video too short: {video_path} "
            f"(needs {end_frame} frames, has {total_frames})"
        )
        cap.release()
        return filename, []

    # Seek to start frame
    cap.set(cv2.CAP_PROP_POS_FRAMES, start_frame)

    # Video properties
    frame_width = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH))
    frame_height = int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))

    # Output video writer (doubled width to show raw + mask side-by-side)
    out = cv2.VideoWriter(
        output_file,
        cv2.VideoWriter_fourcc(*"XVID"),
        fps,
        (frame_width * 2, frame_height),
    )

    bg_subtractor = cv2.createBackgroundSubtractorMOG2(
        history=1000, varThreshold=15, detectShadows=False
    )

    frame_buffer = []  # For temporal frame averaging
    all_prey_counts = []

    current_frame = start_frame

    while current_frame < end_frame:
        ret, frame = cap.read()
        if not ret:
            break

        # Convert to grayscale and blur
        gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
        gray = cv2.GaussianBlur(gray, (5, 5), 1)

        
        # ---------------
        # Prey Mask
        # ---------------

        # Maintain a rolling buffer of grayscale frames
        frame_buffer.append(gray)
        if len(frame_buffer) > FRAME_BUFFER_INTERVAL:
            frame_buffer.pop(0)

        # Average the buffer to reduce temporal noise        
        average_frame = np.mean(frame_buffer, axis=0).astype(np.uint8)
        
        # Apply background subtraction
        prey_mask = bg_subtractor.apply(average_frame)

        # Morphological cleanup: open + close to remove noise
        kernel = np.ones((2, 2), np.uint8)
        prey_mask = cv2.morphologyEx(prey_mask, cv2.MORPH_OPEN, kernel)
        prey_mask = cv2.morphologyEx(prey_mask, cv2.MORPH_CLOSE, kernel)
        
        
        # ---------------
        # Fish Mask
        # ---------------
        
        # Apply adaptative threshold
        fish_mask = cv2.adaptiveThreshold(
            gray,
            255,
            cv2.ADAPTIVE_THRESH_GAUSSIAN_C,
            cv2.THRESH_BINARY,
            101,  # Size of local neighborhood
            10  # Constant subtracted from local mean
        )
        
        # Ensure black background
        fish_mask = cv2.bitwise_not(fish_mask)

        # Morphological cleanup: open + close to remove noise
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (5, 5))
        fish_mask = cv2.morphologyEx(fish_mask, cv2.MORPH_CLOSE, kernel)
        fish_mask = cv2.morphologyEx(fish_mask, cv2.MORPH_OPEN, kernel)
        
        # Remove small objects
        num_labels, labels, stats, _ = cv2.connectedComponentsWithStats(
            fish_mask, connectivity=8
        )

        fish_mask_cleaned = np.zeros_like(fish_mask)

        for i in range(1, num_labels):  # skip background (label 0)
            area = stats[i, cv2.CC_STAT_AREA]
            if area >= MAX_PREY_SIZE:
                fish_mask_cleaned[labels == i] = 255
                
        # Dilate AFTER cleaning
        kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (7, 7))
        fish_mask_cleaned = cv2.morphologyEx(fish_mask_cleaned, cv2.MORPH_DILATE, kernel)

        
        # -----------------
        # Merged Mask
        # -----------------
        
        # Merge the prey mask and the fish mask
        merged_mask = cv2.bitwise_or(prey_mask, fish_mask_cleaned)
        
        # Remove small holes
        kernel = np.ones((2, 2), np.uint8)
        merged_mask = cv2.morphologyEx(merged_mask, cv2.MORPH_CLOSE, kernel)
        
        contours, _ = cv2.findContours(
            merged_mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE
        )

        prey_count = 0
        merged_colored = cv2.cvtColor(merged_mask, cv2.COLOR_GRAY2BGR)  # For visualization

        for contour in contours:
            area = cv2.contourArea(contour)
            if MIN_PREY_SIZE < area < MAX_PREY_SIZE:
                prey_count += 1

                # Get centroid of the contour
                M = cv2.moments(contour)
                if M["m00"] != 0:
                    cx = int(M["m10"] / M["m00"])
                    cy = int(M["m01"] / M["m00"])

                    # Draw count on both raw and binary images
                    cv2.putText(
                        frame, str(prey_count), (cx, cy),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1
                    )
                    cv2.putText(
                        merged_colored, str(prey_count), (cx, cy),
                        cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1
                    )

        # Display current frame's prey count
        cv2.putText(
            frame, f"Prey Count: {prey_count}", (10, 15),
            cv2.FONT_HERSHEY_SIMPLEX, 0.5, (0, 255, 255), 1
        )

        # Update prey count history
        all_prey_counts.append(prey_count)
        
        # Combine frame with binary mask visualization and save
        combined_frame = np.hstack((frame, merged_colored))
        out.write(combined_frame)

        current_frame += 1

    cap.release()
    out.release()
    
    if len(all_prey_counts) != window_frames:
        print(f"Warning: {filename} returned {len(all_prey_counts)} frames")

    return filename, all_prey_counts

# --- File discovery ---
files = [f for f in os.listdir(input_path) 
         if "diff" not in f and f not in existing_filenames]

if len(files) == 0:
    print("No video found.")
else:
    print("Videos found:")
    for f in files:
        print(f)
        
    # --- Extract metadata from file names ---
    date_arr = []
    video_num_arr = []
    type_arr = []
    lux_arr = []
    dpf_arr = []
    richness_arr = []
    delay_arr = []
    
    lux_pattern = r"\d+(?:\.\d+)?lux"
    dpf_pattern = r"\d+dpf"
    delay_pattern = r"\d+start"
    
    for file_name in files:
        date_arr.append(file_name[:10])
        video_num_arr.append(int(file_name[11:13]))
    
        if "MeCP2" in file_name[14:19]:
            type_arr.append("MeCP2")
        elif "WT" in file_name[14:19]:
            type_arr.append("WT")
        else:
            type_arr.append(np.nan)
    
        lux_m = re.search(lux_pattern, file_name)
        lux_arr.append(float(lux_m.group()[:-3]) if lux_m else np.nan)
    
        dpf_m = re.search(dpf_pattern, file_name)
        dpf_arr.append(int(dpf_m.group()[:-3]) if dpf_m else np.nan)
    
        if "peu" in file_name:
            richness_arr.append("low")
        elif "bcp" in file_name:
            richness_arr.append("high")
        else:
            richness_arr.append(np.nan)
            
        delay_m = re.search(delay_pattern, file_name)
        delay_arr.append(int(delay_m.group()[:-5]) if delay_m else None)
    
    video_df = pd.DataFrame({
        'file_name': files,
        'date': date_arr,
        'video_num': video_num_arr,
        'type': type_arr,
        'lux': lux_arr,
        'dpf': dpf_arr,
        'richness': richness_arr,
    })
    
    video_df = pd.concat([existing_video_df, video_df], ignore_index=True)
    video_df.to_csv(video_data_csv, index=False)
    print("Video metadata appended and saved.")
    
    # --- Prepare for multiprocessing ---
    args = [(os.path.join(input_path, file), output_path, delay_arr[i]) 
            for i, file in enumerate(files)]
    
    # Process videos in parallel, collect results
    with Pool(processes=os.cpu_count()) as pool:
        results = list(
            tqdm(pool.imap_unordered(process_video_wrapper, args), 
                 total=len(args), desc="Processing Videos")
        )
    print("All videos processed.")
    
    count_dict = {}
    
    for filename, all_prey_counts in results:
        if not all_prey_counts:
            continue  # skipped / too short / failed
    
        count_dict[filename] = all_prey_counts
        
    lengths = {k: len(v) for k, v in count_dict.items()}
    if len(set(lengths.values())) != 1:
        raise ValueError(f"Frame count mismatch: {lengths}")
    
    count_df = pd.DataFrame.from_dict(count_dict, orient="columns")
    if not existing_count_df.empty:
        count_df = pd.concat([existing_count_df, count_df], axis=1)
    count_df = count_df.loc[:, ~count_df.columns.duplicated()]
    count_df.to_csv(counts_csv, index=False)
    print("Prey counts appended and saved.")