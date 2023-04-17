#!/bin/bash



# define array of launch files
#file_names=("10" "20" "30" "40" "50" "60" "70" "80" "90")
#file_names=("10" "20" "25" "33" "50")

# loop over each method and downsample percent
#for name in "${file_names[@]}"
#do
    # make new folder
    new_folder="exact"
    echo $new_folder
    mkdir "$new_folder"

    exec_name="gpr_mapper_bag"
    echo $exec_name

    # time the process, run the ros executable
    start=$(date +%s.N)
    sleep 5
    rosrun sensor_stream_ros $exec_name "wiggles_curve_main.bag" "$new_folder"
    end=$(date +%s.%N)
    duration=$(echo "$end - $start" | bc)
    echo $duration > "$new_folder.txt"
done

