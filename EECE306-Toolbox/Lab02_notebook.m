%% EECE 306 Lab 2 notebook. Coordinate Systems, Points and Vector Components
% *Team 8.* Thalea Collymore, Kroix Jones, Xavier Moore
%
% Fill in every TODO, run |publish('Lab02_notebook.m')| from the toolbox
% root, print the HTML to PDF, three to five pages.

%% Setup
clear; close all;
rng(306);

%% Point roundtrips, 1000 random points
r = randn(1000,3);
esph = max(em.vec.mag(em.coord.sph2c(em.coord.c2sph(r)) - r));
ecyl = max(em.vec.mag(em.coord.cyl2c(em.coord.c2cyl(r)) - r));
fprintf('spherical   roundtrip max error = %.3e\n', esph);
fprintf('cylindrical roundtrip max error = %.3e\n', ecyl);

%% Poles and the phi convention
disp('c2sph([0 0 1])  ='); disp(em.coord.c2sph([0 0 1]))
disp('c2sph([0 0 -1]) ='); disp(em.coord.c2sph([0 0 -1]))
q = [ 1 1 0; -1 1 0; -1 -1 0; 1 -1 0 ];
p = em.coord.c2cyl(q);
fprintf('phi in the four quadrants = %.4f %.4f %.4f %.4f rad\n', p(:,2));
fprintf('all inside [0, 2*pi) = %d\n', all(p(:,2) >= 0 & p(:,2) < 2*pi));

%% Vector component roundtrips and the invariant
A  = randn(1000,3);
As = em.coord.vecC2Sph(A, r);
ev = max(em.vec.mag(em.coord.vecSph2C(As, r) - A));
em_ = max(abs(em.vec.mag(As) - em.vec.mag(A)));
fprintf('vector roundtrip (spherical pair) max error = %.3e\n', ev);
fprintf('magnitude change under conversion max       = %.3e\n', em_);
Ac = em.coord.vecC2Cyl(A, r);
ev2 = max(em.vec.mag(em.coord.vecCyl2C(Ac, r) - A));
fprintf('vector roundtrip (cylindrical pair) max error = %.3e\n', ev2);

%% The two position demonstration
% The same Cartesian vector A = 2x + 3y, expressed in spherical
% components at two different positions.
A0 = [2 3 0];
disp('at (0, 5, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [0 5 0]))
disp('at (4, 0, 0), (Ar, Atheta, Aphi) ='); disp(em.coord.vecC2Sph(A0, [4 0 0]))
% The vector stayed at 2x +3x. The shperical coordinates changed because the spherical unit vectors are dependent on their position
%A position arguement is mandatory for component converson as it determines the oreintation of the spherical coordinates at the location. 
% Point conversion is meaningless as they are the same and not relative to the coordinate directions.

%% Interpretation
% the roundtrip tests catches a bug where the two functions are inconsistant with each other, which can happen from  incorrect formulas,
%or coordinates not in the right space. It does not catch if both functions make the same mistake, as converting back and forth would 
%would still give the same point. The invariant that catches this error is the vector magnitude. It is refernce free because
%because it does not depend on any coordiante systems. 

%% Problems encountered
% Some vectors modules didnt run correctly the first time due to either incorrect variables used or typos.

%% Full test suite
runTests
