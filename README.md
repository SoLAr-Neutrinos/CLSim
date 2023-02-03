# CLSim
Parametrized charge and light simulation

## Install 

1. Make sure root6 is installed and active 
2. Run `make` 

## Input Files 

Use [qpixg4](https://github.com/Q-Pix/qpixg4) for the G4 generation.  (We should probably fork that at some point to have something like this - Talk to Johnny about this!)

## Run 

`./analyze_light InputOfG4 PlacementForPD PDSize OutputFile (--charge PlacementForCD --pixSize CDSize --number NEvents)`

For more details look at the Docs folder


## HEP MC CLUSTER 

As we need to use `screen` as a submissio system, there is a example submit script in the `RunFiles` folder. 
Change to your setup and submit several jobs to the background. 
If to many cores are occupied already, these jobs get crashed by the kernel. 
So watch out what you are doing. Its not a replacement for a proper submission system, just a workaround as such a system is absent... 