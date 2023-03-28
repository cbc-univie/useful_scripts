"""
velreporter.py: Outputs simulation trajectories in DCD format

"""
from __future__ import absolute_import
__author__ = "Peter Eastman"
__version__ = "1.0"

import array
import os
import time
import struct
import math
import simtk.openmm as mm
from simtk.unit import picoseconds, nanometers, angstroms, is_quantity, norm
from simtk.openmm import Vec3
from simtk.openmm.app.internal.unitcell import computeLengthsAndAngles

class VELFile(object):
    """VELFile provides methods for creating DCD VELfiles."""

    def __init__(self, file, topology, dt, firstStep=0, interval=1, append=False):
        """Create a DCD file and write out the header, or open an existing file to append.

        Parameters
        ----------
        file : file
            A file to write to
        topology : Topology
            The Topology defining the molecular system being written
        dt : time
            The time step used in the trajectory
        firstStep : int=0
            The index of the first step in the trajectory
        interval : int=1
            The frequency (measured in time steps) at which states are written
            to the trajectory
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        """
        self._file = file
        self._topology = topology
        self._firstStep = firstStep
        self._interval = interval
        self._modelCount = 0
        if is_quantity(dt):
            dt = dt.value_in_unit(picoseconds)
        dt /= 0.04888821
        self._dt = dt
        boxFlag = 0
        if topology.getUnitCellDimensions() is not None:
            boxFlag = 1
        if append:
            file.seek(8, os.SEEK_SET)
            self._modelCount = struct.unpack('<i', file.read(4))[0]
            file.seek(268, os.SEEK_SET)
            numAtoms = struct.unpack('<i', file.read(4))[0]
            if numAtoms != len(list(topology.atoms())):
                raise ValueError('Cannot append to a DCD file that contains a different number of atoms')
        else:
            header = struct.pack('<i4c9if', 84, b'C', b'O', b'R', b'D', 0, firstStep, interval, 0, 0, 0, 0, 0, 0, dt)
            header += struct.pack('<13i', boxFlag, 0, 0, 0, 0, 0, 0, 0, 0, 24, 84, 164, 2)
            header += struct.pack('<80s', b'Created by OpenMM')
            header += struct.pack('<80s', b'Created '+time.asctime(time.localtime(time.time())).encode('ascii'))
            header += struct.pack('<4i', 164, 4, len(list(topology.atoms())), 4)
            file.write(header)


    def writeModelVel(self, positions, unitCellDimensions=None, periodicBoxVectors=None):
        """Write out a model to the DCD file.

        The periodic box can be specified either by the unit cell dimensions
        (for a rectangular box), or the full set of box vectors (for an
        arbitrary triclinic box).  If neither is specified, the box vectors
        specified in the Topology will be used. Regardless of the value
        specified, no dimensions will be written if the Topology does not
        represent a periodic system.

        Parameters
        ----------
        positions : list
            The list of atomic positions to write
        unitCellDimensions : Vec3=None
            The dimensions of the crystallographic unit cell.
        periodicBoxVectors : tuple of Vec3=None
            The vectors defining the periodic box.
        """
        if len(list(self._topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
        if is_quantity(positions):
            positions = positions.value_in_unit(nanometers/picoseconds)
        if any(math.isnan(norm(pos)) for pos in positions):
            raise ValueError('Particle position is NaN')
        if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError('Particle position is infinite')
        file = self._file

        self._modelCount += 1
        if self._interval > 1 and self._firstStep+self._modelCount*self._interval > 1<<31:
            # This will exceed the range of a 32 bit integer.  To avoid crashing or producing a corrupt file,
            # update the header to say the trajectory consisted of a smaller number of larger steps (so the
            # total trajectory length remains correct).
            self._firstStep //= self._interval
            self._dt *= self._interval
            self._interval = 1
            file.seek(0, os.SEEK_SET)
            file.write(struct.pack('<i4c9if', 84, b'C', b'O', b'R', b'D', 0, self._firstStep, self._interval, 0, 0, 0, 0, 0, 0, self._dt))

        # Update the header.

        file.seek(8, os.SEEK_SET)
        file.write(struct.pack('<i', self._modelCount))
        file.seek(20, os.SEEK_SET)
        file.write(struct.pack('<i', self._firstStep+self._modelCount*self._interval))

        # Write the data.

        file.seek(0, os.SEEK_END)
        boxVectors = self._topology.getPeriodicBoxVectors()
        if boxVectors is not None:
            if periodicBoxVectors is not None:
                boxVectors = periodicBoxVectors
            elif unitCellDimensions is not None:
                if is_quantity(unitCellDimensions):
                    unitCellDimensions = unitCellDimensions.value_in_unit(nanometers)
                boxVectors = (Vec3(unitCellDimensions[0], 0, 0), Vec3(0, unitCellDimensions[1], 0), Vec3(0, 0, unitCellDimensions[2]))*nanometers
            (a_length, b_length, c_length, alpha, beta, gamma) = computeLengthsAndAngles(boxVectors)
            a_length = a_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            b_length = b_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            c_length = c_length * 10.  # computeLengthsAndAngles returns unitless nanometers, but need angstroms here.
            angle1 = math.sin(math.pi/2-gamma)
            angle2 = math.sin(math.pi/2-beta)
            angle3 = math.sin(math.pi/2-alpha)
            file.write(struct.pack('<i6di', 48, a_length, angle1, b_length, angle2, angle3, c_length, 48))
        length = struct.pack('<i', 4*len(positions))
        for i in range(3):
            file.write(length)
            data = array.array('f', (10*x[i] for x in positions))
            data.tofile(file)
            file.write(length)
        try:
            file.flush()
        except AttributeError:
            pass

class VELReporter(object):
    """VELReporter outputs a series of frames from a Simulation to a DCD file.

    To use it, create a VELReporter, then add it to the Simulation's list of reporters.
    """

    def __init__(self, file, reportInterval, append=False, enforcePeriodicBox=None):
        """Create a VELReporter.

        Parameters
        ----------
        file : string
            The file to write to
        reportInterval : int
            The interval (in time steps) at which to write frames
        append : bool=False
            If True, open an existing DCD file to append to.  If False, create a new file.
        enforcePeriodicBox: bool
            Specifies whether particle positions should be translated so the center of every molecule
            lies in the same periodic box.  If None (the default), it will automatically decide whether
            to translate molecules based on whether the system being simulated uses periodic boundary
            conditions.
        """
        self._reportInterval = reportInterval
        self._append = append
        self._enforcePeriodicBox = enforcePeriodicBox
        if append:
            mode = 'r+b'
        else:
            mode = 'wb'
        self._out = open(file, mode)
        self._dcd = None

    def describeNextReport(self, simulation):
        """Get information about the next report this object will generate.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for

        Returns
        -------
        tuple
            A six element tuple. The first element is the number of steps
            until the next report. The next four elements specify whether
            that report will require positions, velocities, forces, and
            energies respectively.  The final element specifies whether
            positions should be wrapped to lie in a single periodic box.
        """
        steps = self._reportInterval - simulation.currentStep%self._reportInterval
        return (steps, False, True, False, False, self._enforcePeriodicBox)

    def report(self, simulation, state):
        """Generate a report.

        Parameters
        ----------
        simulation : Simulation
            The Simulation to generate a report for
        state : State
            The current state of the simulation
        """

        if self._dcd is None:
            self._dcd = VELFile(
                self._out, simulation.topology, simulation.integrator.getStepSize(),
                simulation.currentStep, self._reportInterval, self._append
            )
        self._dcd.writeModelVel(state.getVelocities(), periodicBoxVectors=state.getPeriodicBoxVectors())

    def __del__(self):
        self._out.close()
