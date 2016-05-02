SmartWire branch
=========

SmartWire branch contains an attempt to have one interchangible library to support I2C on chips with USI, full TWI, or soft-I2C only with one library. The goal would be for no code changes to be needed inside libraries in order to support I2C master functionality, because it's an unpleasant barrier to have to modify a library (which sounds scary - though it's not that hard, it's still a pain to have to do) 

I2C slave functionality is not supported - since you don't typically have libraries that depend on using Wire in slave mode, this is not a priority.

This code is not expected to be usable in it's current state
