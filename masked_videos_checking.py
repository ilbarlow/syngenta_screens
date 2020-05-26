#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  5 16:46:34 2020

@author: ibarlow

Script for checking wells:
    1. Get a list of all the masked videos for a single day of imaging

    2. Find unique plates as a single plate is imaged three times to result
    in 3x videos; we only need to mark the bad wells on one of these videos

    3. Launch the tierpsy_gui multiwormtracker to mark the bad wells

    4. Save and then move on to the next video

"""

from pathlib import Path
import pandas as pd
import re
import numpy as np

# day_directories = [Path('/Volumes/SAMSUNG/SyngentaStrainScreen/MaskedVideos/20200306')]
# INPUT_DIR = Path('/Volumes/SAMSUNG/SyngentaStrainScreen/MaskedVideos')
INPUT_DIR = Path('/Volumes/behavgenom$/Ida/Data/Hydra/SyngentaScreen/' +
                'MaskedVideos')

N_CAMERAS = 6
run_regex = r"(?<=run)\d{1,2}"
OVERWRITE = False
imgstore_regex = r"\d{8,}(?=/metadata)"

if __name__ == '__main__':
    day_directories = [f for f in list(INPUT_DIR.glob('*')) if f.is_dir()]

    for day in day_directories:
        aux_day_dir = Path(str(day).replace('MaskedVideos',
                           'AuxiliaryFiles'))
        if not aux_day_dir.exists():
            aux_day_dir.mkdir()

        masked_videos_outfile = aux_day_dir /\
            '{}_masked_videos_well_checking_second_checks.csv'.format(day.stem)

        if not OVERWRITE:
            if masked_videos_outfile.exists():
                print('{} file already exists, nothing done here')

        try:
            masked_videos = list(day.rglob('metadata.hdf5'))
            prestim_videos = [f for f in masked_videos if
                              'prestim' in str(f)]

            # put the prestim videos in a df for exporting
            prestim_videos = pd.DataFrame({'masked_video': prestim_videos})
            prestim_videos['run_number'] = [int(re.search(run_regex,
                                            str(r.masked_video))[0])
                                            for i, r in
                                            prestim_videos.iterrows()]
            prestim_videos['imgstore_camera'] = [
                                        int(re.search(imgstore_regex,
                                                      str(r.masked_video))
                                                        [0]) for i, r in
                                                prestim_videos.iterrows()
                                                ]
            imgstore_cameras = list(np.unique(
                    prestim_videos['imgstore_camera'].to_list()))

            if len(imgstore_cameras) < 6:
                print('Not all cameras were recording on this day')
                N_CAMERAS = len(imgstore_cameras)

            # check all the videos are there
            assert prestim_videos.shape[0] % N_CAMERAS == 0

        except AssertionError:
            print('Videos missing on {}'.format(day))

            missing_videos_outfile = aux_day_dir /\
                '{}_missing_masked_videos_v2.csv'.format(day.stem)

            missing_videos = []
            prestim_videos_grouped = prestim_videos.groupby('run_number')
            runs = list(prestim_videos['run_number'].unique())
            for r in runs:
                if prestim_videos_grouped.get_group(r).shape[0] % 6 != 0:
                    _checking = prestim_videos_grouped.get_group(r)
                    missing_cameras = list(set(_checking['imgstore_camera']
                                                ).symmetric_difference(
                                                    imgstore_cameras))
                    missing_videos.append((r,
                                           [str(m) for m in
                                            missing_cameras]))

            missing_videos = pd.DataFrame(missing_videos,
                                          columns=['run_number',
                                                   'camera'])

            missing_videos.to_csv(missing_videos_outfile, index=False)

        prestim_videos.to_csv(masked_videos_outfile, index=False)
