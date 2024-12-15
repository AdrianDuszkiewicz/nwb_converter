from pathlib import Path

import nwb_converter.nwb_functions as nwb_func
import nwb_converter.import_functions as imp_func


def main():
    datapath = Path('/Users/doctordu/Dropbox (Personal)/Professional/Data/Trisha')
    foldername = 'H7101-241201'
    metaname = 'H7100_metadata.xlsx'

    # Load metadata from the xml file you keep your metadata in
    metadata = imp_func.get_metadata(datapath, metaname, foldername)
    # load start time from Metadata.txt file
    start_time = imp_func.get_start_time(datapath, foldername)
    # Load tracking, epochs and spikes from Matlab files
    pos, ang, epochs, spikes, shank_id, waveforms, maxIx, tr2pk = imp_func.get_data_matlab(datapath, foldername)
    # Load event times and process a bit
    events = imp_func.get_events(datapath, foldername, time_offset=epochs.at[2, 'Start'])
    events = events.iloc[16:]  # remove first TTLs when squares switch on. TODO: remove hard-coding

    # create and fill the NWB file
    nwbfile = nwb_func.create_nwb_file(metadata, start_time)
    nwbfile = nwb_func.add_epochs(nwbfile, epochs, metadata)
    nwbfile = nwb_func.add_tracking(nwbfile, pos, ang)
    nwbfile = nwb_func.add_probe(nwbfile, metadata, n_shanks=2, n_channels=32, step=25, electrodes_so_far=0)
    nwbfile = nwb_func.add_units(nwbfile, spikes, waveforms, shank_id)
    nwbfile = nwb_func.add_events(nwbfile, events)

    # save NWB file
    nwb_func.save_nwb_file(nwbfile, datapath, foldername)

    # load NWB file (for testing)
    data = nwb_func.load_nwb_file(datapath, foldername)

    print()

if __name__ == '__main__':
    main()

