# Cellpose workflow

The Cellpose workflow was run with Cellpose 3.1.0 in Python 3.11. The environment specification is provided in `cellpose-environment.yml` and `cellpose-environment-frozen.yml`.

The script processes raw TIFF images exported from Leica LAS X. It processes one emission channel at a time; the channel is selected by the `file_identifier` setting in `run_cellpose.py`. Configure the input location, channel identifier, and Cellpose settings for the dataset being analysed before running the script.
