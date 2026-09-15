classdef femmsession < fpproc
%FEMMSESSION Stateful, in-memory magnetic analysis session.
%   S = FEMMSESSION(FILENAME) loads a magnetic .fem problem and owns its
%   native model, mesher, solver, trial solution, and accepted state. Native
%   pointers never leave session_interface_mex; OBJECTHANDLE is an opaque
%   class handle used only by this wrapper.
%
%   A typical evaluation is:
%       s = xfemm.femmsession('motor.fem');
%       s.setBackend('tangle');
%       s.setCircuit('phase-a', 'current', 10);
%       s.setAGEPosition('airgap', rotorAngle, 0);
%       status = s.solve();
%       s.accept();                 % or s.reject()
%
%   The accepted state can optionally be persisted with SAVESTATE and
%   restored into another session for the same model with LOADSTATE.
%   After SOLVE, the object is also a fully initialized FPPROC instance;
%   all magnetic point, contour, block, circuit, mesh, and air-gap analysis
%   methods documented by FPPROC can be called directly on the session. The
%   post-processor is refreshed in memory; SOLVE does not write an .ans file.
%
%   FEMMSESSION methods:
%     setBackend      - select 'triangle' (default) or 'tangle' meshing
%     mesh             - create the mesh now and return its element count
%     setCircuit       - set current/voltage/open/coupled circuit constraint
%     setAGEPosition   - set inner and outer AGE angles in degrees
%     setFrequency     - set analysis frequency
%     setTime          - set the evaluation time
%     setRotationalInstances - mesh one tile and repeat it by rotation
%     instancedInfo    - identities, instance count, materialisation count
%     setInstanceOverride - per-instance circuit/turns/magnetisation override
%     solve            - solve and return compact status/statistics
%     result           - explicitly return the complete latest trial data
%     accept           - accept the trial as the next initial state
%     reject           - discard the trial without advancing state
%     saveState        - save the accepted state to a MAT-file
%     loadState        - restore an accepted state from a MAT-file

    properties (Access = private, Hidden = true)
        sessionHandle
    end

    methods
        function this = femmsession(filename)
            %FEMMSESSION Load FILENAME and create a native analysis session.
            narginchk(1, 1);
            this@fpproc();
            this.sessionHandle = session_interface_mex('new', filename);
            this.openfilename = filename;
            this.FemmProblem = loadfemmfile(filename);
        end

        function delete(this)
            %DELETE Release every native component owned by the gateway.
            if ~isempty(this.sessionHandle)
                session_interface_mex('delete', this.sessionHandle);
                this.sessionHandle = [];
            end
        end

        function setBackend(this, name)
            %SETBACKEND Select the 'triangle' or 'tangle' meshing backend.
            session_interface_mex('backend', this.sessionHandle, lower(name));
        end

        function count = mesh(this)
            %MESH Create/recreate the session mesh and return element count.
            count = session_interface_mex('mesh', this.sessionHandle);
        end

        function setCircuit(this, name, constraint, value)
            %SETCIRCUIT Set a named circuit constraint for the next solve.
            %   Constraint is 'current', 'voltage', 'open', or 'coupled'. A
            %   real or complex VALUE is required for current and voltage.
            if nargin < 4, value = 0; end
            session_interface_mex('circuit', this.sessionHandle, name, ...
                                   lower(constraint), value);
        end

        function setAGEPosition(this, name, innerAngle, outerAngle)
            %SETAGEPOSITION Set named AGE inner/outer angles, in degrees.
            session_interface_mex('age', this.sessionHandle, name, ...
                                   innerAngle, outerAngle);
        end

        function setFrequency(this, frequency)
            %SETFREQUENCY Set frequency in hertz for subsequent solves.
            session_interface_mex('frequency', this.sessionHandle, frequency);
        end

        function setTime(this, time)
            %SETTIME Set the user-defined time associated with the trial.
            session_interface_mex('time', this.sessionHandle, time);
        end

        function setRotationalInstances(this, centerX, centerY, instanceCount, totalAngle)
            %SETROTATIONALINSTANCES Mesh one tile and repeat it by rotation.
            %   SETROTATIONALINSTANCES(CX, CY, N, ANGLE) meshes the loaded
            %   problem once as a tile and places N instances rotated about
            %   (CX, CY) so that they cover ANGLE degrees (default 360). The
            %   tile's matched periodic boundaries become the welded seams.
            if nargin < 5, totalAngle = 360; end
            session_interface_mex('instance', this.sessionHandle, centerX, ...
                                   centerY, instanceCount, totalAngle);
        end

        function info = instancedInfo(this)
            %INSTANCEDINFO Return instancing identities and counters.
            info = session_interface_mex('instanceinfo', this.sessionHandle);
        end

        function setInstanceOverride(this, instance, sourceLabel, circuit, ...
                                     magnetisationRotation, currentScale)
            %SETINSTANCEOVERRIDE Set per-instance physics for one block label.
            %   SETINSTANCEOVERRIDE(INSTANCE, SOURCELABEL, CIRCUIT, MAGROT, SCALE)
            %   uses 1-based indices. CIRCUIT 0 leaves circuit membership
            %   unchanged; NaN MAGROT/SCALE leave those quantities unchanged.
            if nargin < 6, currentScale = NaN; end
            if nargin < 5, magnetisationRotation = NaN; end
            if nargin < 4, circuit = 0; end
            session_interface_mex('instanceoverride', this.sessionHandle, ...
                                   instance - 1, sourceLabel - 1, circuit - 1, ...
                                   magnetisationRotation, currentScale);
        end

        function out = solve(this)
            %SOLVE Synchronize, solve, and return compact status information.
            %   The returned struct contains success, solution identity and
            %   time, mesh sizes, and solver execution counters. Field and
            %   circuit data are deliberately not copied into this result;
            %   request quantities with the inherited FPPROC methods instead.
            out = session_interface_mex('solve', this.sessionHandle);
            this.isdocopen = true;
        end

        function out = result(this)
            %RESULT Explicitly return the complete latest trial data.
            %   This includes nodal A and coordinates and can be large. Prefer
            %   the inherited FPPROC analysis methods for requested quantities.
            out = session_interface_mex('result', this.sessionHandle);
        end

        function accept(this)
            %ACCEPT Make the latest trial the initial state for future solves.
            session_interface_mex('accept', this.sessionHandle);
        end

        function reject(this)
            %REJECT Discard the latest trial and retain the accepted state.
            session_interface_mex('reject', this.sessionHandle);
        end

        function saveState(this, filename)
            %SAVESTATE Persist the current accepted state in a MAT-file.
            state = session_interface_mex('state', this.sessionHandle); %#ok<NASGU>
            if isempty(state), error('MFEMM:session:noState', 'No solution has been accepted.'); end
            save(filename, 'state');
        end

        function loadState(this, filename)
            %LOADSTATE Restore an accepted state previously saved by SAVESTATE.
            saved = load(filename, 'state');
            if ~isfield(saved, 'state')
                error('MFEMM:session:badStateFile', 'File does not contain a session state.');
            end
            session_interface_mex('restore', this.sessionHandle, saved.state);
        end
    end

    methods (Access = protected)
        function varargout = callPostProcessor(this, varargin)
            [varargout{1:nargout}] = session_interface_mex(varargin{1}, ...
                                                            this.sessionHandle, ...
                                                            varargin{2:end});
        end
    end
end
