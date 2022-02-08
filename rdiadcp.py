import numpy as np


def adcp_ascii(filnam):
    """
    % ADCP_ASCII read classic WinRiver TXT ASCII files
    % input: filename (string)
    % output: ADCP structure
    % file structure is based on:
    % "WinRiver II User's Guide" (Teledyne RD Instruments, Feb. 2007) pp. 64-65
    %
    % no error checking is done, so this could fail miserably on malformed
    % files
    %
    % D. Nowacki nowacki@uw.edu Feb. 2012
    %
    % Revisions:
    % 28 May 2014 - more comments
    % 7 April 2021 - convert to Python -- D. Nowacki dnowacki@usgs.gov.
    """
    adcp = {}
    with open(filnam) as f:

        # these three lines are at the beginning of each file.
        adcp["comment1"] = f.readline()
        adcp["comment2"] = f.readline()
        line = f.readline().split()

        (
            adcp["binsize"],
            adcp["blank"],
            adcp["firstbin"],
            adcp["numcells"],
            adcp["numpings"],
            adcp["tpe"],
            adcp["watermode"],
        ) = [int(x) for x in line]

        thevars = [
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
            "hunsec",
            "ensnum",
            "numens",
            "pitch",
            "roll",
            "heading",
            "temp",
            "btvel_x",
            "btvel_y",
            "btvel_z",
            "btvel_err",
            "depthsounder",
            "ggaalt",
            "ggadeltaalt",
            "ggahdop",
            "depth1",
            "depth2",
            "depth3",
            "depth4",
            "elapdist",
            "elaptime",
            "distnorth",
            "disteast",
            "distmadegood",
            "lat",
            "lon",
            "qmid",
            "qtop",
            "qbot",
            "binstofollow",
            "measunit",
            "velref",
            "intensunits",
            "intensscale",
            "soundabsorp",
            "z",
            "spd",
            "dir",
            "east",
            "north",
            "up",
            "err",
            "echo1",
            "echo2",
            "echo3",
            "echo4",
            "percgood",
            "q",
        ]
        for var in thevars:
            adcp[var] = []

        # now read blocks throughout file
        # note that, to conserve space, ensembles are indexed sequentially, not by
        # ensemble number.

        line = f.readline().split()
        ens = 0
        while line:
            # ROW 1
            # R Field Description
            # 1 1     ENSEMBLE TIME -Year (at start of ensemble)
            #   2                   - Month
            #   3                   - Day
            #   4                   - Hour
            #   5                   - Minute
            #   6                   - Second
            #   7                   - Hundredths of seconds
            #   8     ENSEMBLE NUMBER (or SEGMENT NUMBER for processed or averaged raw  data)
            #   9     NUMBER OF ENSEMBLES IN SEGMENT (if averaging ON or processing data)
            #   10    PITCH ñ Average for this ensemble (degrees)
            #   11    ROLL ñ Average for this ensemble (degrees)
            #   12    CORRECTED HEADING - Average ADCP heading (corrected for one cycle error) + heading offset + magnetic variation
            #   13    ADCP TEMPERATURE - Average for this ensemble (∞C)

            thevars = [
                "year",
                "month",
                "day",
                "hour",
                "minute",
                "second",
                "hunsec",
                "ensnum",
                "numens",
                "pitch",
                "roll",
                "heading",
                "temp",
            ]
            for n, var in enumerate(thevars):
                adcp[var].append(float(line[n]))

            # ROW 2
            # 2	1   BOTTOM-TRACK VELOCITY - East(+)/West(-); average for this ensemble (cm/s or ft/s)
            #   2   Reference = BTM       - North(+)/South(-)
            #   3                         - Vertical (up[+]/down[-])
            #   4                         - Error
            # 2 1   BOTTOM-TRACK VELOCITY ñ GPS (GGA or VTG) Velocity (calculated from GGA String) Reference = GGA  East(+)/West (-1)
            #   2   Reference = VTG - GPS (GGA or VTG) North(+)/South(-) Velocity
            #   3                         - BT  (up[+]/down[-]) Velocity
            #   4                         - BT Error

            line = f.readline().split()
            thevars = [
                "btvel_x",
                "btvel_y",
                "btvel_z",
                "btvel_err",
                "depthsounder",
                "ggaalt",
                "ggadeltaalt",
                "ggahdop",
                "depth1",
                "depth2",
                "depth3",
                "depth4",
            ]
            for n, var in enumerate(thevars):
                adcp[var].append(float(line[n]))

            # ROW 3

            line = f.readline().split()
            thevars = ["elapdist", "elaptime", "distnorth", "disteast", "distmadegood"]
            for n, var in enumerate(thevars):
                adcp[var].append(float(line[n]))

            # ROW 4

            line = f.readline().split()
            thevars = [
                "lat",
                "lon",
            ]  # dont care about other values on this line for now
            for n, var in enumerate(thevars):
                adcp[var].append(float(line[n]))

            # ROW 5

            line = f.readline().split()
            thevars = [
                "qmid",
                "qtop",
                "qbot",
            ]  # dont care about other values on this line for now
            for n, var in enumerate(thevars):
                adcp[var].append(float(line[n]))

            # ROW 6
            # 6 1	NUMBER OF BINS TO FOLLOW
            #   2   MEASUREMENT UNIT ñ cm or ft
            #   3   VELOCITY REFERENCE ñ BT, GGA, VTG, or NONE for current velocity data rows 7-26 fields 2-7
            #   4   INTENSITY UNITS - dB or counts
            #   5   INTENSITY SCALE FACTOR ñ in dB/count
            #   6   SOUND ABSORPTION FACTOR ñ in dB/m

            line = f.readline().split()
            thevars = [
                "binstofollow",
                "measunit",
                "velref",
                "intensunits",
                "intensscale",
                "soundabsorp",
            ]
            for n, var in enumerate(thevars):
                try:
                    adcp[var].append(float(line[n]))
                except ValueError:
                    adcp[var].append(line[n])

            # ROW 7-26
            # 7-26  1    DEPTH ñ Corresponds to depth of data for present bin (depth cell); includes ADCP depth and blanking value; in m or ft.
            #       2   VELOCITY MAGNITUDE
            #       3   VELOCITY DIRECTION
            #       4   EAST VELOCITY COMPONENT ñ East(+)/West(-)
            #       5   NORTH VELOCITY COMPONENT - North(+)/South(-)
            #       6   VERTICAL VELOCITY COMPONENT - Up(+)/Down(-)
            #       7   ERROR VELOCITY
            #       8   BACKSCATTER     ñ Beam 1
            #       9                   - Beam 2
            #       10                  - Beam 3
            #       11                  - Beam 4
            #       12  PERCENT-GOOD
            #       13  DISCHARGE

            thevars = [
                "z",
                "spd",
                "dir",
                "east",
                "north",
                "up",
                "err",
                "echo1",
                "echo2",
                "echo3",
                "echo4",
                "percgood",
                "q",
            ]
            tmp = {}
            for var in thevars:
                tmp[var] = []
            for n in range(adcp["numcells"]):
                line = f.readline().split()
                for n, var in enumerate(thevars):
                    tmp[var].append(float(line[n]))
            for n, var in enumerate(thevars):
                adcp[var].append(tmp[var])

            line = f.readline().split()
            ens += 1
            # break

    # do brief QAQC on the data
    for k in adcp:
        if isinstance(adcp[k], list):
            adcp[k] = np.array(adcp[k])
    for k in ["spd", "east", "north", "up", "err", "lat", "lon"]:
        adcp[k][adcp[k] == -32768] = np.nan
    for k in ["lat", "lon"]:
        adcp[k][adcp[k] == 30000] = np.nan

    return adcp
