function [x, y] = radial_machine_sample_points ()
%RADIAL_MACHINE_SAMPLE_POINTS Deterministic field sample points for comparison.
%
%   [x, y] = radial_machine_sample_points ()
%
% Returns points in the two meshed halves of the air gap of the checked-in
% 12-pole/36-slot design (the region between the two AGE arcs is bridged by the
% air-gap element and is deliberately not meshed). The seed is fixed so the
% redraw, sliding, and tiled cases sample identical locations.

    Rmo = 0.0565;
    g = 2e-3;
    rng (20240919);
    n = 8;
    margin = 0.1 * g;
    lowerRadii = Rmo + margin + (g/3 - 2*margin) * rand (n, 1);
    upperRadii = Rmo + 2*g/3 + margin + (g/3 - 2*margin) * rand (n, 1);
    radii = [lowerRadii; upperRadii];
    angles = 60 * rand (2*n, 1);
    x = radii .* cosd (angles);
    y = radii .* sind (angles);
end
