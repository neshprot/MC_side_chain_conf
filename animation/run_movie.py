import os
import re
import time


if __name__ == '__main__':

    def atoi(text):
        return int(text) if text.isdigit() else text


    def natural_keys(text):
        return [atoi(c) for c in re.split('(\d+)', text)]


    def next_frame(i, frame):
        os.system(f'write_tcl.sh {i} {frame}')


    i = 0
    frames = []
    directory = 'frames'
    for filename in os.scandir(directory):
        if filename.is_file():
            frames.append(filename.path[7::])
    frames.sort(key=natural_keys)

    for i, frame in enumerate(frames):
        next_frame(i, frame)
