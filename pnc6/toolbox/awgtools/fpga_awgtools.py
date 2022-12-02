"""
Provides experiment classes that add AWG use to the FPGA.
"""
import time
import logging

import numpy as np

import awgloader
import pulseseq.sequencer as sequencer
import pulseseq.pulselib as pulselib

from fpga_lib.parameters import *
from fpga_lib.dsl import *


class AWGExperimentAddon(Parameterized):
    awg = None
    awg_clock = 1e9
    do_setup_awg = True
    awg_ssb_calibrations = {}
    awg_timeout = 10
    load_awg_sequence = BoolParameter(True)


    def __init__(self, *arg, **kw):
        super(AWGExperimentAddon, self).__init__(*arg, **kw)


    def awg_sequence(self):
        raise NotImplementedError


    def awg_file_loader(self):
        raise NotImplementedError


    def setup_awg(self, run=True, wait=True):
        timeout = self.awg_timeout
        awg = self.instruments[self.awg]

        if self.load_awg_sequence:
            awg.set_clock(self.awg_clock)
            s = self.awg_sequence()
            s = sequencer.Sequencer(s)

            ssb = []
            for _caln in self.awg_ssb_calibrations:
                ssb.append(sequencer.SSB(if_period=self.awg_ssb_calibrations[_caln]['if_period'],
                                         chans=self.awg_ssb_calibrations[_caln]['chans'],
                                         dphi=self.awg_ssb_calibrations[_caln]['dphi'],
                                         amps=self.awg_ssb_calibrations[_caln]['amps'],
                                         outchans=self.awg_ssb_calibrations[_caln]['outchans']))
            
            for _ssb in ssb:
                s.add_ssb(_ssb)

            for i in range(1,5):
                s.add_required_channel(i)
            
            seqs = s.render()
            l = self.awg_file_loader()
            l.load(seqs)
            
            sequencer.Pulse.clear_pulse_data()
        
        if run:
            awg.run()
            if wait:
                if hasattr(self, 'logger'):
                    self.logger.info("Waiting for AWG... (timeout = {} s)".format(timeout))
                else:
                    logging.info("Waiting for AWG... (timeout = {} s)".format(timeout))
                
                if timeout > 0:
                    t0 = time.time()
                    waitingtime = 0
                    i = 0
                    while (time.time() < t0 + timeout):
                        waitingtime = time.time() - t0
                        if waitingtime > i*10:
                            print "Waiting... ({} s)".format(i*10)
                            i += 1
                        try:
                            rs = awg.get_runstate()
                        except:
                            rs = False
                        if rs:
                            break
                        time.sleep(0.05)
                    if not rs:
                        raise Exception("AWG not starting. Aborting experiment.")
                    else:
                        logging.info("AWG seems ready... wait another second to make sure.")
                        time.sleep(1)

# class AWGExperimentAddon