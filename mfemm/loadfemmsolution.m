function [FemmProblem, Solution] = loadfemmsolution(filename, problemonly)
% loads a femm solution from a file and stores it in a FemmProblem and
% FemmProblem solution structure
%
% Syntax
%
% [FemmProblem, Solution] = loadfemmsolution(filename)
% [FemmProblem, Solution] = loadfemmsolution(filename, problemonly)
%
% Input
%
%   filename - string containing the location of the femm solution file to
%     be loaded
%
%   problemonly - flag to determine whether the solution is actually
%     loaded, or only the FemmProblem description. Defaults to false,
%     meaning everything is loaded. If true, Solution will be empty.
%
% Output
%
%   FemmProblem - a structure describing the finite element problem, see
%     writefemmfile.m for details.
%
%   Solution - sturucture containing the following fields:
%
%       Info - Structure containing information about the solution. It
%         contains three fields: NNodes, NElements and NCircuits, which are
%         the number of nodes, triangle elements and circuits in the
%         solution respectively.
%
%       MeshNodes - Matrix of node information, where each row represents a
%         node. If the problem frequency is zero, it is an (n x 3) matrix
%         where the first column is the x (or r) coordinate of each node,
%         and the second the y (or z) coordinate. The thrd column in this
%         case is the value of the vector potential (A) calculated at this
%         node. For axisymmetric problems, this is actually 2 x pi x r x A,
%         as the vector potential is calculated on a contour at the radius
%         of the node.
%
%         If the problem frequency is not zero, this will be an (n x 4)
%         matrix where the final two columns are the real and imaginary
%         components of the vector potential.
%
%       Vertices - (n x 3) matrix of mesh vertices. Each row list the
%         indices of nodes at the three vertices of the mesh triangle
%         element.
%
%       LabelNum - (n x 1) vector of label numbers associated with each
%         element.
%
%       Circuits - array of structures containing information on each
%         circuit found in the solution. Each structure has the following
%         fields: Type, Volts, CurrentDensity
%
%

% Copyright 2012 Richard Crozier
% 
%    Licensed under the Apache License, Version 2.0 (the "License");
%    you may not use this file except in compliance with the License.
%    You may obtain a copy of the License at
% 
%        http:%www.apache.org/licenses/LICENSE-2.0
% 
%    Unless required by applicable law or agreed to in writing, software
%    distributed under the License is distributed on an "AS IS" BASIS,
%    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
%    See the License for the specific language governing permissions and
%    limitations under the License.

    if nargin < 2
        problemonly = false;
    end
    
    % use a structure to create an enum of file types
    ftype.magnetics = 0;
    ftype.heatflow = 1;

    % use a structure to map the file extensions to problem types
    fexmap.ans = 'magnetics';
    fexmap.anh = 'heatflow';
    fexmap.fem = 'magnetics';
    fexmap.feh = 'heatflow';
    
    % try to determine from the file name extension
    [~, ~, fext] = fileparts (filename);


    if ~isempty(fext) && isfield(fexmap, fext(2:end))
        FemmProblem.ProbInfo.Domain = fexmap.(fext(2:end));
    else
        error ('MFEMM:badfiletype', ...
               'The output file name extension is not recognised and you have not explicitly set the file type.');
    end
    
%     FemmProblem.ProbInfo.PrevSolutionFile = '';
%     FemmProblem.ProbInfo.PrevSolutionType = 0;
    
    fid = fopen(filename);

    cleanupObj = onCleanup(@() fclose(fid));

    matind = 0;
    boundind = 0;

    q = fgetl(fid);

    while ischar(q)

        % trim whitespace from the string at start and end
        q = strtrim(q);

        % Deal with flag for file format version
        if strncmpi(q,'[format]',8)

            v = str2double(StripKey(q));

            switch ftype.(FemmProblem.ProbInfo.Domain)
                
                case ftype.magnetics
                    if v ~= 4
                        error('This file is from a different version of FEMM\nRe-analyze the problem using the current version.');
                    end
                    
                case ftype.heatflow
                    if v ~= 1
                        error('This file is from a different version of FEMM\nRe-analyze the problem using the current version.');
                    end
                    
                otherwise
                    error('Unrecognised file extension.');
                    
            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % Frequency of the problem
        if strncmpi(q,'[frequency]',11)
            v = StripKey(q);
            FemmProblem.ProbInfo.Frequency = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        % Depth of the problem
        if strncmpi(q,'[depth]',7)
            v = StripKey(q);
            FemmProblem.ProbInfo.Depth = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        % Minimum angle constraint
        if strncmpi(q,'[MinAngle]',7)
            v = StripKey(q);
            FemmProblem.ProbInfo.MinAngle = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end
        
        % Option to force use of default max mesh, overriding
        % user choice
        if strncmpi(q,'[forcemaxmesh]',14)
            v = StripKey(q);
            % 0 == do not override user mesh choice
            % not 0 == do override user mesh choice
            if str2double(v) == 0
                FemmProblem.ProbInfo.ForceMaxMeshArea = false;
            else
                FemmProblem.ProbInfo.ForceMaxMeshArea = true;
            end
            q = fgetl(fid);
            continue;
        end

        % Option to use smart meshing
        if strncmpi(q,'[dosmartmesh]',13)
            v = StripKey(q);
            % 0 == do not use smart mesh
            % not 0 == use smart mesh
            if str2double(v) == 0
                FemmProblem.ProbInfo.SmartMesh = false;
            else
                FemmProblem.ProbInfo.SmartMesh = true;
            end
            q = fgetl(fid);
            continue;
        end
        
        % Units of length used by the problem
        if strncmpi(q,'[lengthunits]',13)

            v = StripKey(q);

            if strncmpi(v,'inches',6)
                FemmProblem.ProbInfo.LengthUnits = 0;
            elseif strncmpi(v,'millimeters',11)
                FemmProblem.ProbInfo.LengthUnits = 1;
            elseif strncmpi(v,'centimeters',1)
                FemmProblem.ProbInfo.LengthUnits = 2;
            elseif strncmpi(v,'mils',4)
                FemmProblem.ProbInfo.LengthUnits  = 4;
            elseif strncmpi(v,'microns',6)
                FemmProblem.ProbInfo.LengthUnits  = 5;
            elseif strncmpi(v,'meters',6)
                FemmProblem.ProbInfo.LengthUnits  = 3;
            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end
        
        if strncmpi(q,'[precision]',11)
            
            v = StripKey(q);
            
            FemmProblem.ProbInfo.Precision = str2double(v);
            
            % get the next line of input
            q = fgetl(fid);
            continue;
            
        end

        % Problem Type (planar or axisymmetric)
        if strncmpi(q,'[problemtype]',13)

            v = StripKey(q);

            if strncmpi(v,'planar',6)
                FemmProblem.ProbInfo.ProblemType = 0;
            end

            if strncmpi(v,'axisymmetric',3)
                FemmProblem.ProbInfo.ProblemType = 1;
            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % Coordinates (cartesian or polar)
        if strncmpi(q,'[coordinates]',13)

            v = StripKey(q);

            if strncmpi(v,'cartesian',4)
                FemmProblem.ProbInfo.Coords = 0;
            end

            if strncmpi(v,'polar',5)
                FemmProblem.ProbInfo.Coords = 1;
            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end
        
        % Choice of solver
        if strncmpi(q,'[ACSolver]',10)

            v = StripKey(q);

            FemmProblem.ProbInfo.ACSolver = str2double(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end
        
        % Previous solution file
        if strncmpi(q,'[prevsoln]',8)

            v = StripKey(q);

            FemmProblem.ProbInfo.PrevSolutionFile = removequotes(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end
        
        % previous solution type
        if strncmpi(q,'[prevtype]',8)

            v = StripKey(q);

            FemmProblem.ProbInfo.PrevSolutionType = str2double(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % Comments
        if (strncmpi(q,'[comment]',9))

            v = StripKey(q);

            FemmProblem.ProbInfo.Comments = removequotes(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % properties for axisymmetric external region
        if( strncmpi(q,'[extzo]',7))
            v = StripKey(q);
            FemmProblem.ProbInfo.extZo = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'[extro]',7))
            v = StripKey(q);
            FemmProblem.ProbInfo.extRo = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'[extri]',7))
            v = StripKey(q);
            FemmProblem.ProbInfo.extRi = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        % Point Properties
        if( strncmpi(q,'<beginpoint>',11))

            % Make a new empty point property
            [FemmProblem, ppropind] = addpointprop_mfemm (FemmProblem, 'New Point Property');

            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<pointname>',11))

            v = StripKey(q);

            FemmProblem.PointProps(ppropind).Name = removequotes(v);

            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<A_re>',6))
            v = StripKey(q);
            FemmProblem.PointProps(ppropind).A_re = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<A_im>',6))
            v = StripKey(q);
            FemmProblem.PointProps(ppropind).A_im = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<I_re>',6))
            v = StripKey(q);
            FemmProblem.PointProps(ppropind).I_re = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<I_im>',6))
            v = StripKey(q);
            FemmProblem.PointProps(ppropind).I_im = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<endpoint>',9))
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        % Boundary Properties;
        if( strncmpi(q,'<beginbdry>',11))

            [FemmProblem, boundind] = addboundaryprop_mfemm(FemmProblem, 'New Boundary', 0);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<bdryname>',10))

            v = StripKey(q);

            FemmProblem.BoundaryProps(boundind).Name = removequotes(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<bdrytype>',10))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).BdryType = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<mu_ssd>',8))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).Mu_ssd = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<sigma_ssd>',11))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).Sigma_ssd = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<A_0>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).A0 = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<A_1>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).A1 = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<A_2>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).A2 = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<phi>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).Phi = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<c0>',4))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).c0 = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<c1>',4))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).c1 = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<c0i>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).c0i = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<c1i>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).c1i = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end
        
        if( strncmpi(q,'<innerangle>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).InnerAngle = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<outerangle>',5))
            v = StripKey(q);
            FemmProblem.BoundaryProps(boundind).OuterAngle = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<endbdry>',9))
            % get the next line of input
            q = fgetl(fid);
            continue;
        end


        % Block Properties;
        if( strncmpi(q,'<beginblock>',12))

            [FemmProblem, matind] = addmaterials_mfemm(FemmProblem, emptymatstruct_mfemm());
            
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<blockname>',10))
            v = StripKey(q);

            FemmProblem.Materials(matind).Name = removequotes(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<mu_x>',6))
            v = StripKey(q);
            FemmProblem.Materials(matind).Mu_x = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<mu_y>',6))
            v = StripKey(q);
            FemmProblem.Materials(matind).Mu_y = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<H_c>',5))
            v = StripKey(q);
            FemmProblem.Materials(matind).H_c = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<J_re>',6))
            v = StripKey(q);
            FemmProblem.Materials(matind).J_re = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<J_im>',6))

            v = StripKey(q);

            if (FemmProblem.ProbInfo.Frequency ~= 0)
                FemmProblem.Materials(matind).J_im = str2double(v);
            end
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<sigma>',7))

            v = StripKey(q);
            FemmProblem.Materials(matind).Sigma = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<phi_h>',7))

            v = StripKey(q);
            FemmProblem.Materials(matind).Phi_h = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<phi_hx>',8))

            v = StripKey(q);
            FemmProblem.Materials(matind).Phi_hx = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<phi_hy>',8))

            v = StripKey(q);
            FemmProblem.Materials(matind).Phi_hy = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<d_lam>',7))

            v = StripKey(q);
            FemmProblem.Materials(matind).d_lam = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<LamFill>',8))

            v = StripKey(q);
            FemmProblem.Materials(matind).LamFill = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<LamType>',9))

            v = StripKey(q);
            FemmProblem.Materials(matind).LamType = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<NStrands>',10))

            v = StripKey(q);
            FemmProblem.Materials(matind).NStrands = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<WireD>',7))

            v = StripKey(q);
            FemmProblem.Materials(matind).WireD = str2double(v);
            % get the next line of input
            q = fgetl(fid);
            continue;
        end

        if( strncmpi(q,'<BHPoints>',10))

            v = StripKey(q);

            nbhpoints = str2double(v);

            if nbhpoints > 0

                for j = 1:nbhpoints

                    s = fgetl(fid);

                    FemmProblem.Materials(matind).BHPoints(j,:) = sscanf(s,'%f');

                end

            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<endblock>',9))

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % Circuit Properties
        if( strncmpi(q,'<begincircuit>',14))

            [FemmProblem, circind] = addcircuit_mfemm(FemmProblem, '');

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<circuitname>',13))

            v = StripKey(q);

            FemmProblem.Circuits(circind).Name = removequotes(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<totalamps_re>',14))

            v = StripKey(q);

            FemmProblem.Circuits(circind).TotalAmps_re = str2double(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<totalamps_im>',14))

            v = StripKey(q);

            if FemmProblem.ProbInfo.Frequency ~= 0
                FemmProblem.Circuits(circind).TotalAmps_im = str2double(v);
            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<circuittype>',13))

            v = StripKey(q);

            FemmProblem.Circuits(circind).CircType = str2double(v);

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if( strncmpi(q,'<endcircuit>',12))

            % get the next line of input
            q = fgetl(fid);
            continue;

        end
        
%         if strncmpi(q,'[BdryProps]',10)
%             
%             v = StripKey(q);
% 
%             nbdrys = round(str2double(v));
%             
%             for i = 1:nbdrys
% 
%                 fprintf(fp,'  <BeginBdry>\n');
%                 fprintf(fp,'    <BdryName> = \"%s\"\n', FemmProblem.BoundaryProps(i).Name);
%                 fprintf(fp,'    <BdryType> = %i\n',FemmProblem.BoundaryProps(i).BdryType);
%                 fprintf(fp,'    <A_0> = %.17g\n',FemmProblem.BoundaryProps(i).A_0);
%                 fprintf(fp,'    <A_1> = %.17g\n',FemmProblem.BoundaryProps(i).A_1);
%                 fprintf(fp,'    <A_2> = %.17g\n',FemmProblem.BoundaryProps(i).A_2);
%                 fprintf(fp,'    <Phi> = %.17g\n',FemmProblem.BoundaryProps(i).Phi);
%                 fprintf(fp,'    <c0> = %.17g\n',FemmProblem.BoundaryProps(i).c0);
%                 fprintf(fp,'    <c0i> = %.17g\n',FemmProblem.BoundaryProps(i).c0i);
%                 fprintf(fp,'    <c1> = %.17g\n',FemmProblem.BoundaryProps(i).c1);
%                 fprintf(fp,'    <c1i> = %.17g\n',FemmProblem.BoundaryProps(i).c1i);
%                 fprintf(fp,'    <Mu_ssd> = %.17g\n',FemmProblem.BoundaryProps(i).Mu_ssd);
%                 fprintf(fp,'    <Sigma_ssd> = %.17g\n',FemmProblem.BoundaryProps(i).Sigma_ssd);
%                 fprintf(fp,'  <EndBdry>\n');
% 
%             end
% 
%         end
        
        % Points list;
        if(strncmpi(q,'[numpoints]',11))

            v = StripKey(q);

            npoints = round(str2double(v));

            for i = 1:npoints

                s = fgetl(fid);
                
                switch ftype.(FemmProblem.ProbInfo.Domain)
                    case ftype.magnetics
                        C = textscan(s, '%f %f %f %f');
                    case ftype.heatflow
                        C = textscan(s, '%f %f %f %f %f');
                        nodeprops.InConductor = C{5};
                end
                
                nodeprops.InGroup = C{4};

                if C{3} > 0

                    nodeprops.PointPropName = FemmProblem.PointProps(C{3}).Name;

                end

                FemmProblem = addnodes_mfemm(FemmProblem, C{1}, C{2}, nodeprops);

                clear nodeprops;

            end
            
            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % read in segment list
        if(strncmpi(q,'[numsegments]',13))

            v = StripKey(q);

            nsegments = round(str2double(v));

            for i = 1:nsegments

                s = fgetl(fid);

                switch ftype.(FemmProblem.ProbInfo.Domain)
                    case ftype.magnetics
                        C = textscan(s, '%f %f %f %f %f %f');
                    case ftype.heatflow
                        C = textscan(s, '%f %f %f %f %f %f %f');
                        segprops.InConductor = C{7};
                end

                segprops.MaxSideLength = C{3};

                segprops.BoundaryMarker = C{4};
                
                if segprops.BoundaryMarker == 0
                    segprops.BoundaryMarker = '';
                else
                    segprops.BoundaryMarker = FemmProblem.BoundaryProps(segprops.BoundaryMarker).Name;
                end
                
                segprops.Hidden = C{5};
                
                segprops.InGroup = C{6};

                FemmProblem = addsegments_mfemm(FemmProblem, C{1}, C{2}, segprops);

                clear segprops;

            end
            
            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % read in arc segment list
        if(strncmpi(q,'[numarcsegments]',16))
            v = StripKey(q);

            narcsegments = round(str2double(v));

            for i = 1:narcsegments

                s = fgetl(fid);

                % fixed: moved here to be set before MeshedSideLength default load!
                arcsegprops.MaxSegDegrees = C{4};
                
                switch ftype.(FemmProblem.ProbInfo.Domain)
                    
                    case ftype.magnetics
                        
                        C = textscan(s,'%f %f %f %f %f %f %f %f'); 
                        
                        if ~isempty (C{8})
                            % old femm versions didn't have the 8th parameter
                            if C{8} > 0
                                arcsegprops.MeshedSideLength = C{8};
                            else
                                arcsegprops.MeshedSideLength = arcsegprops.MaxSegDegrees; % ###note: same as MeshedSideLength by default... 
                            end
                        end

                    case ftype.heatflow
                        
                        C = textscan(s,'%f %f %f %f %f %f %f %f');
                        
                        arcsegprops.InConductor = C{8};
                end
                
                

                arcsegprops.BoundaryMarker = C{5};
                
                % zero boundary marker corresponds to no boundary marker
                if arcsegprops.BoundaryMarker == 0
                    arcsegprops.BoundaryMarker = '';
                else
                    arcsegprops.BoundaryMarker = FemmProblem.BoundaryProps(arcsegprops.BoundaryMarker).Name;
                end

                arcsegprops.Hidden = C{6};
                
                arcsegprops.InGroup = C{7};
                
                FemmProblem = addarcsegments_mfemm(FemmProblem, C{1}, C{2}, C{3}, arcsegprops);

                clear arcsegprops;

            end

            % get the next line of input
            q = fgetl(fid);
            continue;
        end


        % read in list of holes;
        if(strncmpi(q,'[numholes]',10))

            v = StripKey(q);

            nholes = round(str2double(v));

            if nholes > 0

                for i = 1:nholes

                    s = fgetl(fid);

                    C = textscan(s,'%f %f %d');

                    FemmProblem.Holes(i) = struct ( 'Coords', [C{1}, C{2}], ...
                                                    'InGroup', C{3} );

                    %  don't add holes to the list
                    %  of block labels because it messes up the
                    %  number of block labels.
                end
            else
                FemmProblem.Holes = [];
            end
            
            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        % read in regional attributes
        if(strncmpi(q,'[numblocklabels]',16))

            v = StripKey(q);

            nregionalattrib = str2double(v);

            for i = 1:nregionalattrib

                s = fgetl(fid);
                
                switch ftype.(FemmProblem.ProbInfo.Domain)
                    
                    case ftype.magnetics
                        
                        C = textscan(s,'%f %f %f %f %f %f %f %f %f %s');

                        % copy the data to the new block
                        FemmProblem.BlockLabels(i).Coords = [C{1}, C{2}];
                        FemmProblem.BlockLabels(i).BlockType = C{3};
                        FemmProblem.BlockLabels(i).MaxArea = C{4};
                        FemmProblem.BlockLabels(i).InCircuit = C{5};
                        FemmProblem.BlockLabels(i).MagDir = C{6};
                        FemmProblem.BlockLabels(i).InGroup = C{7};
                        FemmProblem.BlockLabels(i).Turns = C{8};
                        FemmProblem.BlockLabels(i).IsExternal = C{9};
                        FemmProblem.BlockLabels(i).IsDefault = bitand (int32 (FemmProblem.BlockLabels(i).IsExternal), int32 (2));
                        FemmProblem.BlockLabels(i).IsExternal = bitand (int32 (FemmProblem.BlockLabels(i).IsExternal), int32 (1));
                        
                        if ~isempty (C{10})
                            % textscan keeps the surrounding quotes, but
                            % writefemmfile adds its own; strip them so a
                            % loadfemmfile/writefemmfile round-trip does not
                            % turn "theta" into ""theta"".
                            fct = strtrim (C{10}{1});
                            if numel (fct) >= 2 && fct(1) == '"' && fct(end) == '"'
                                fct = fct(2:end-1);
                            end
                            FemmProblem.BlockLabels(i).MagDirFctn = fct;
                        else
                            FemmProblem.BlockLabels(i).MagDirFctn = '';
                        end

                        if FemmProblem.BlockLabels(i).MaxArea < 0
                            FemmProblem.BlockLabels(i).MaxArea = 0;
        %                 else
        %                     % convert from a side length to area constraint
        %                     FemmProblem.BlockLabels(i).MaxArea = pi * FemmProblem.BlockLabels(i).MaxArea * FemmProblem.BlockLabels(i).MaxArea / 4;
                        end
                        
                        % ###fix: this one is expected in addblocklabel_mfemm(), but not set anywhere when problem loaded from file (only set when problem generated from scratch)!
                        FemmProblem.BlockLabels(i).IsMagnet = ~isnan(FemmProblem.BlockLabels(i).MagDir); 
                
                    case ftype.heatflow
                        
                        C = textscan(s,'%f %f %f %f %f %f');
                        
                        FemmProblem.BlockLabels(i).Coords = [C{1}, C{2}];
                        FemmProblem.BlockLabels(i).BlockType = C{3};
                        FemmProblem.BlockLabels(i).MaxArea = C{4};
                        FemmProblem.BlockLabels(i).InGroup = C{5};
                        FemmProblem.BlockLabels(i).IsDefault = C{6};
                        
                end

            end

            % get the next line of input
            q = fgetl(fid);
            continue;

        end

        if strncmpi(q,'[solution]',10) 
            hassolution = true;
            break;
        else
            hassolution = false;
            % get the next line of input
            q = fgetl(fid);
            continue;
        end
        


    end
    
%     if isfield(FemmProblem, 'Segments')
%         % convert the boudary marker number in the segment list to names
%         for i = 1:numel(FemmProblem.Segments)
%             
%             if FemmProblem.Segments(i).BoundaryMarker == 0
%                 FemmProblem.Segments(i).BoundaryMarker = '';
%             else
%                 FemmProblem.Segments(i).BoundaryMarker = ...
%                     FemmProblem.BoundaryProps(FemmProblem.Segments(i).BoundaryMarker).Name;
%             end
%             
%         end
%     end
%     
%     if isfield(FemmProblem, 'ArcSegments')
%         % convert the boudary marker number in the arc segment list to names
%         for i = 1:numel(FemmProblem.ArcSegments)
%             
%             if FemmProblem.ArcSegments(i).BoundaryMarker == 0
%                 FemmProblem.ArcSegments(i).BoundaryMarker = '';
%             else
%                 FemmProblem.ArcSegments(i).BoundaryMarker = ...
%                     FemmProblem.BoundaryProps(FemmProblem.ArcSegments(i).BoundaryMarker).Name;
%             end
%             
%         end
%     end

    if isfield(FemmProblem, 'BlockLabels')
        for i = 1:numel(FemmProblem.BlockLabels)
            FemmProblem.BlockLabels(i).BlockType = FemmProblem.Materials(FemmProblem.BlockLabels(i).BlockType).Name;
        end
    end
    
    switch ftype.(FemmProblem.ProbInfo.Domain)
        
        case ftype.magnetics
            
            if isfield(FemmProblem, 'BlockLabels')
                for i = 1:numel(FemmProblem.BlockLabels)
                    if FemmProblem.BlockLabels(i).InCircuit > 0
                        % convert the circuit number to a circuit name
                        FemmProblem.BlockLabels(i).InCircuit = FemmProblem.Circuits(FemmProblem.BlockLabels(i).InCircuit).Name;
                    end
                end
            end
    end

    if hassolution && ~problemonly

        % read in mesh node locations and nodal values

        % get the number of nodes to read in
        s = fgetl(fid);
        
        Solution.Info.NNodes = cell2mat(textscan(s,'%f\n'));
        
        for i = 1:Solution.Info.NNodes

            % get the string containing the node info
            s = fgetl(fid);

            if FemmProblem.ProbInfo.Frequency ~= 0
                % read in four values from the string, the x and y
                % coordinates of the node and the real and imaginary
                % values of A at the node
                Solution.MeshNodes(i,1:4) = sscanf(s,'%f %f %f %f');
            else
                % read in three values from the string, the x and y
                % coordinates of the node and the value of A at the node
                Solution.MeshNodes(i,1:3) = sscanf(s,'%f %f %f');
                Solution.MeshNodes(i,4) = 0;
            end

        end

        % read in the element data

        % get the number of elements to read in
        s = fgetl(fid);
        Solution.Info.NElements = cell2mat(textscan(s,'%f\n'));

        for i = 1:Solution.Info.NElements

            s = fgetl(fid);

            vals = sscanf(s,'%f %f %f %f');

            Solution.Vertices(i,:) = vals(1:3);
            Solution.LabelNum(i,1) = vals(4);
%             Solution.LabelStr{i} = FemmProblem.BlockLabels(vals(4)+1).BlockType;

        end

        % get the number of circuits to read in
        s = fgetl(fid);
        Solution.Info.NCircuits = cell2mat(textscan(s,'%f'));

        for i = 1:Solution.Info.NCircuits
            
            if i == 1
                Solution.Circuits = struct('Type', [], 'Volts', [], 'CurrentDensity', []);
            end

            s = fgetl(fid);

            if FemmProblem.ProbInfo.Frequency == 0

                cvals = sscanf(s,'%f %f');

                Solution.Circuits(i).Type = cvals(1);

                if Solution.Circuits(i).Type == 0
                    Solution.Circuits(i).Volts = cvals(2);
                else
                    Solution.Circuits(i).CurrentDensity = cvals(2);
                end

            else

                cvals = sscanf(s,'%f %f %f');

                Solution.Circuits(i).Type = cvals(1);

                if Solution.Circuits(i).Type == 0
                    Solution.Circuits(i).Volts = cvals(2) + 1j*cvals(3);
                else
                    Solution.Circuits(i).CurrentDensity = cvals(2) + 1j*cvals(3);
                end

            end

        end

%         Solution = postprocsolution(FemmProblem, Solution);
        
    else
        Solution = [];
    end
    

end


function v = StripKey(q)

    [start_idx, end_idx, extents, matches, tokens, names, splits] = regexp(q, '(?<=\=).*');

    v = strtrim(matches{1});

end


function str = removequotes(v)

    [tokenstring] = regexp(v, '"(.*)"', 'tokens');

    str = char(tokenstring{1});

end


% function [Solution] = postprocsolution(FemmProblem, Solution)
% 
%     % convert depth from whatever the supplied problem units were to meters 
%     % for internal computations
%     if (FemmProblem.ProbInfo.Depth == -1)
%         FemmProblem.ProbInfo.Depth = 1; 
%     else
%         FemmProblem.ProbInfo.Depth = length2metres_mfemm(FemmProblem.ProbInfo.Depth, FemmProblem.ProbInfo.LengthUnits);
%     end
% 
% 	% Calculate triangle element centroids 
%     % triangle centroid given by = (1/3 (x1+x2+x3), 1/3 (y1+y2+y3))  
%     % see wikipedia page for further info:
%     % http://en.wikipedia.org/wiki/Centroid#Of_triangle_and_tetrahedron
%     for i = 1:size(Solution.Vertices,1)
%         Solution.TriCentroids(i,1:2) = [sum(Solution.MeshNodes(Solution.Vertices(i,:)+1,1)),  ...
%                                         sum(Solution.MeshNodes(Solution.Vertices(i,:)+1,2))];
%     end
%     Solution.TriCentroids = Solution.TriCentroids ./ 3;
%     
%     % Compute the element radii
%     rsqrx = realpow(  reshape(Solution.MeshNodes(Solution.Vertices(:)+1,1), size(Solution.Vertices)) ...
%                     - reshape(Solution.TriCentroids(Solution.Vertices(:)+1,1), size(Solution.Vertices)), ...
%                     2);
%                     
%     rsqry = realpow(  reshape(Solution.MeshNodes(Solution.Vertices(:)+1,2), size(Solution.Vertices)) ...
%                     - reshape(Solution.TriCentroids(Solution.Vertices(:)+1,2), size(Solution.Vertices)), ...
%                     2);                    
%     
%     Solution.TriRadii =  max(rsqrx + rsqry, [], 2);
% 
%     % calculate the magnetization direction for each element 
%     Solution = calcelementmagdir(FemmProblem, Solution);
%    
%     
% % 	% Find flux density in each element;
% % 	for(i=0;i<meshelem.GetSize();i++) GetElementB(meshelem[i]);
% % 
% % 	% Find extreme values of A;
% % 	A_Low=meshnode[0].A.re; A_High=meshnode[0].A.re;
% % 	for(i=1;i<meshnode.GetSize();i++)
% % 	{
% % 		if (meshnode[i].A.re>A_High) A_High=meshnode[i].A.re;
% % 		if (meshnode[i].A.re<A_Low)  A_Low =meshnode[i].A.re;
% % 
% % 		if(Frequency!=0)
% % 		{
% % 			if (meshnode[i].A.im<A_Low)  A_Low =meshnode[i].A.im;
% % 			if (meshnode[i].A.im>A_High) A_High=meshnode[i].A.im;
% % 		}
% % 	}
% % 	% save default values for extremes of A
% % 	A_lb=A_Low;
% % 	A_ub=A_High;
% 
% % 	if(Frequency!=0){ % compute frequency-dependent permeabilities for linear blocks;
% % 
% % 		CComplex deg45; deg45=1+I;
% % 		CComplex K,halflag;
% % 		double ds;
% % 		double w=2.*PI*Frequency;
% % 
% % 		for(k=0;k<blockproplist.GetSize();k++){
% % 		if (blockproplist[k].LamType==0){
% % 			blockproplist[k].mu_fdx = blockproplist[k].mu_x*
% % 									  exp(-I*blockproplist[k].Theta_hx*PI/180.);
% % 			blockproplist[k].mu_fdy = blockproplist[k].mu_y*
% % 									  exp(-I*blockproplist[k].Theta_hy*PI/180.);
% % 
% % 			if(blockproplist[k].Lam_d!=0){
% % 				halflag=exp(-I*blockproplist[k].Theta_hx*PI/360.);
% % 				ds=sqrt(2./(0.4*PI*w*blockproplist[k].Cduct*blockproplist[k].mu_x));
% % 				K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
% % 				if (blockproplist[k].Cduct!=0)
% % 				{
% % 					blockproplist[k].mu_fdx=(blockproplist[k].mu_fdx*tanh(K)/K)*
% % 					blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
% % 				}
% % 				else{
% % 					blockproplist[k].mu_fdx=(blockproplist[k].mu_fdx)*
% % 					blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
% % 				}
% % 
% % 				halflag=exp(-I*blockproplist[k].Theta_hy*PI/360.);
% % 				ds=sqrt(2./(0.4*PI*w*blockproplist[k].Cduct*blockproplist[k].mu_y));
% % 				K=halflag*deg45*blockproplist[k].Lam_d*0.001/(2.*ds);
% % 				if (blockproplist[k].Cduct!=0)
% % 				{
% % 					blockproplist[k].mu_fdy=(blockproplist[k].mu_fdy*tanh(K)/K)*
% % 					blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
% % 				}
% % 				else{
% % 					blockproplist[k].mu_fdy=(blockproplist[k].mu_fdy)*
% % 					blockproplist[k].LamFill+(1.-blockproplist[k].LamFill);
% % 				}
% % 			}
% % 		}
%         
% 
% end
% 
%                 
% function Solution = calcelementmagdir(FemmProblem, Solution)
% % calculate the magnetisation direction in each element
%     
%     Solution.MagDir = zeros(size(Solution.TriRadii, 1), 1);
%     
%     for i = 1:numel(FemmProblem.BlockLabels)
%         if isempty(FemmProblem.BlockLabels(i).MagDirFctn)
%             Solution.MagDir(Solution.LabelNum == i,1) = FemmProblem.BlockLabels(i).MagDir;
%         else
%             x = Solution.TriCentroids(Solution.LabelNum == i,1);
%             y = Solution.TriCentroids(Solution.LabelNum == i,2);
%             % calculate theta in degrees
%             theta = atan2(y,x) .* 180./pi;
%             % calculate the distance from the 
%             R = mfemmdeps.magn([x,y],2);
%             % convert the string containing the magnetization direction
%             % function to an inline matlab function
%             magdirfcn = inline(FemmProblem.BlockLabels(i).MagDirFctn, 'x', 'y', 'theta', 'R');
%             % vectorise the mag direction function
%             magdirfcn = vectorize(magdirfcn);
%             % evaluate the inline function
%             Solution.MagDir(Solution.LabelNum == i,1) = magdirfcn(x, y, theta, R);
%         end
%     end
%     
% end