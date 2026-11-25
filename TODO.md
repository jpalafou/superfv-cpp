[ ] `clip_last_dt` parameter: boolean
[ ] `fprint` parameter: frequency at which a step summary is printed to the command line
        `fprint=0`: print nothing
        `fprint=1`: print every step
        `fprint=n>1` print every n steps, but always print last step
[ ] `nout` parameter: number of outputs to write to disk. requires code for writing conservative cell averages `u` to disk
        `nout=0`: no output
        `nout=1`: write output from last step only
        `nout=2`: write output from initial condition and last step
        `nout>2`: write output from linspace of steps from `t=0` to `t=tstop`
