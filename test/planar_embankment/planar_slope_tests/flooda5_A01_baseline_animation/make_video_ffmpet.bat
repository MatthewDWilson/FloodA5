rem ffmpeg -r 12 -f image2 -i "frame_%%06d.png" -qscale 0 "flooda5_A01_baseline_animation.avi"
ffmpeg -framerate 12 -i "frame_%%06d.png" -vf "scale=1920:1080" -c:v libx264 -crf 20 -pix_fmt yuv420p flooda5_A01_baseline_animation.mp4

