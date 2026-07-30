classdef mfemmsession < handle
%MFEMMSESSION Stateful, in-memory magnetic analysis session.
%   S = MFEMMSESSION(FILENAME) loads a magnetic .fem problem and owns its
%   native model, mesher, solver, trial solution, and accepted state. Native
%   pointers never leave session_interface_mex; OBJECTHANDLE is an opaque
%   class handle used only by this wrapper.
%
%   A typical evaluation is:
%       s = mfemmsession('motor.fem');
%       s.setBackend('tangle');
%       s.setCircuit('phase-a', 'current', 10);
%       s.setAGEPosition('airgap', rotorAngle, 0);
%       result = s.solve();
%       s.accept();                 % or s.reject()
%
%   The accepted state can optionally be persisted with SAVESTATE and
%   restored into another session for the same model with LOADSTATE.
%
%   MFEMMSESSION methods:
%     setBackend      - select 'triangle' (default) or 'tangle' meshing
%     mesh             - create the mesh now and return its element count
%     setCircuit       - set current/voltage/open/coupled circuit constraint
%     setAGEPosition   - set inner and outer AGE angles in degrees
%     setFrequency     - set analysis frequency
%     setTime          - set the evaluation time
%     solve            - create and return a disposable trial result
%     result           - return the most recent trial result
%     accept           - accept the trial as the next initial state
%     reject           - discard the trial without advancing state
%     saveState        - save the accepted state to a MAT-file
%     loadState        - restore an accepted state from a MAT-file

    properties (Access = private, Hidden = true)
        objectHandle
    end

    methods
        function this = mfemmsession(filename)
            %MFEMMSESSION Load FILENAME and create a native analysis session.
            narginchk(1, 1);
            this.objectHandle = session_interface_mex('new', filename);
        end

        function delete(this)
            %DELETE Release every native component owned by the gateway.
            if ~isempty(this.objectHandle)
                session_interface_mex('delete', this.objectHandle);
                this.objectHandle = [];
            end
        end

        function setBackend(this, name)
            %SETBACKEND Select the 'triangle' or 'tangle' meshing backend.
            session_interface_mex('backend', this.objectHandle, lower(name));
        end

        function count = mesh(this)
            %MESH Create/recreate the session mesh and return element count.
            count = session_interface_mex('mesh', this.objectHandle);
        end

        function setCircuit(this, name, constraint, value)
            %SETCIRCUIT Set a named circuit constraint for the next solve.
            %   Constraint is 'current', 'voltage', 'open', or 'coupled'. A
            %   real or complex VALUE is required for current and voltage.
            if nargin < 4, value = 0; end
            session_interface_mex('circuit', this.objectHandle, name, ...
                                   lower(constraint), value);
        end

        function setAGEPosition(this, name, innerAngle, outerAngle)
            %SETAGEPOSITION Set named AGE inner/outer angles, in degrees.
            session_interface_mex('age', this.objectHandle, name, ...
                                   innerAngle, outerAngle);
        end

        function setFrequency(this, frequency)
            %SETFREQUENCY Set frequency in hertz for subsequent solves.
            session_interface_mex('frequency', this.objectHandle, frequency);
        end

        function setTime(this, time)
            %SETTIME Set the user-defined time associated with the trial.
            session_interface_mex('time', this.objectHandle, time);
        end

        function out = solve(this)
            %SOLVE Synchronize, solve, and return a trial result struct.
            %   Fields include circuit current, flux linkage, optional
            %   terminal voltage, nodal magnetic vector potential A, and
            %   corresponding mesh-node coordinates x and y.
            out = session_interface_mex('solve', this.objectHandle);
        end

        function out = result(this)
            %RESULT Return the latest trial without solving again.
            out = session_interface_mex('result', this.objectHandle);
        end

        function accept(this)
            %ACCEPT Make the latest trial the initial state for future solves.
            session_interface_mex('accept', this.objectHandle);
        end

        function reject(this)
            %REJECT Discard the latest trial and retain the accepted state.
            session_interface_mex('reject', this.objectHandle);
        end

        function saveState(this, filename)
            %SAVESTATE Persist the current accepted state in a MAT-file.
            state = session_interface_mex('state', this.objectHandle); %#ok<NASGU>
            if isempty(state), error('MFEMM:session:noState', 'No solution has been accepted.'); end
            save(filename, 'state');
        end

        function loadState(this, filename)
            %LOADSTATE Restore an accepted state previously saved by SAVESTATE.
            saved = load(filename, 'state');
            if ~isfield(saved, 'state')
                error('MFEMM:session:badStateFile', 'File does not contain a session state.');
            end
            session_interface_mex('restore', this.objectHandle, saved.state);
        end
    end
end
