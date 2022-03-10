function PBMtoDEM(PSD, I_max, rad_range)
PSD_frac = PSD/sum(PSD);
N_p = round(sum(PSD));


fileID = fopen('in.txt','w');

fprintf(fileID,'%12s\n','atom_style 	granular');
fprintf(fileID,'%12s\n','atom_modify	map array');
fprintf(fileID,'%12s\n','boundary	f f f');
fprintf(fileID,'%12s\n','newton		off');

fprintf(fileID,'%12s\n','communicate	single vel yes');

fprintf(fileID,'%12s\n','region		reg block -0.0001 0.200001 -0.0001 0.1001 -0.0001 0.1001 units box');
fprintf(fileID,'%12s\n','create_box	1 reg');

fprintf(fileID,'%12s\n','neighbor	0.002 bin');
fprintf(fileID,'%12s\n','neigh_modify	delay 0');

%%%Material properties required for new pair styles

fprintf(fileID,'%12s\n','fix 		m1 all property/global youngsModulus peratomtype 5.e6');
fprintf(fileID,'%12s\n','fix 		m2 all property/global poissonsRatio peratomtype 0.25');
fprintf(fileID,'%12s\n','fix 		m3 all property/global coefficientRestitution peratomtypepair 1 0.1');
fprintf(fileID,'%12s\n','fix 		m4 all property/global coefficientFriction peratomtypepair 1 0.5');
fprintf(fileID,'%12s\n','fix        id all property/global coefficientRollingFriction peratomtypepair 1 0.01');

%%%%%%New pair style%%%%%%%

fprintf(fileID,'%12s\n','pair_style gran model hertz tangential history rolling_friction cdt');
fprintf(fileID,'%12s\n','pair_coeff	* *');

fprintf(fileID,'%12s\n','timestep	0.00001');

fprintf(fileID,'%12s\n','fix		gravi all gravity 9.81 vector 0.0 0.0 -1.0');

%%%%%%%%%%%%%%%%%%%granular walls%%%%%%%%%%%%%%%%%%%%%
fprintf(fileID,'%12s\n','fix		drum all mesh/surface file meshes/drum1.stl type 1 scale 0.2');
fprintf(fileID,'%12s\n','fix 		wall all wall/gran model hertz tangential history rolling_friction cdt mesh n_meshes 1 meshes drum');

%%%%%%%%%%%%%%%%%%%distributions  insertion%%%%%%%%%%%%%%%%
prime=[10007,10009,10037,10039,10061,10067,10069,10079,10091,10093,...
       10099,10103,10111,10133,10139,10141,10151,10159,10163,10169,...
       10177,10181,10193,10211,10223,10243,10247,10253,10259,10267,...
       10271,10273,10289,10301,10303,10313,10321,10331,10333,10337];
for i=1:(I_max)
   % prime = nextprime(sym(i*10000)); % prime number
   fprintf(fileID,'%12s\n',['fix  pts',num2str(i),' all particletemplate/sphere ',...
    num2str(prime(i)),' atom_type 1 density constant 1000 radius constant ',num2str(rad_range(i))]);
end

fprintf(fileID,['fix		pdd1 all particledistribution/discrete/numberbased 32452843 ',num2str(I_max)]);
for i=1:I_max
    fprintf(fileID,[' pts',num2str(i),' ',num2str(PSD_frac(i))]);
end
%%%%%%%%%%%%%%%%%%%region  insertion%%%%%%%%%%%%%%
fprintf(fileID,'\n\n');
fprintf(fileID,'%12s\n','region 		BC1 cylinder x 0.05 0.05 0.05 0.0 0.2 units box');
fprintf(fileID,'%12s\n','group		nve_group region BC1');

%%%%%%%%%%%%%%%%%%deprecated pour command%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% fix		ins nve_group pour/dev mass 30. 1 distributiontemplate pdd1 vol 0.25 200 massflowrate 30. vel uniform 0. 0. 0. 0. 0.0 region bc

%%%%%%%%%%%%%%%%%%%%%%5particle insertion%%%%%%%%%%%
fprintf(fileID,'%12s\n',['fix	ins nve_group insert/pack seed 67867979 distributiontemplate pdd1 insert_every once overlapcheck yes all_in yes vel constant 0. 0. 0. particles_in_region ',num2str(N_p),' region BC1']);


%%%%%%%%%%%%%%%%%%%apply nve integration to all particles that are inserted as single particles%%%%%%%
fprintf(fileID,'%12s\n','fix		integr nve_group nve/sphere');

% 

%%%%%%%%%%%%%%%%%%%%output settings, include total thermal energy%%%%%%%%%%%%%%%%%

fprintf(fileID,'%12s\n','fix				ts all check/timestep/gran 1000 0.1 0.1');
fprintf(fileID,'%12s\n','compute 	contactsearch all contact/atom');
fprintf(fileID,'%12s\n','thermo_style	custom step atoms ke'); %  c_rke f_ts[1] f_ts[2] vol
fprintf(fileID,'%12s\n','thermo			100');
fprintf(fileID,'%12s\n','thermo_modify	lost ignore norm no');

%%%%%%%%%%%%%%%%insert the first particles so that dump is not empty%%%%%%%%%%%%%%%%%
fprintf(fileID,'%12s\n','run			1');

fprintf(fileID,'%12s\n','dump		dmp all custom/vtk 1000 post/conveyor.vtk id type type x y z ix iy iz vx vy vz fx fy fz omegax omegay omegaz radius ');
fprintf(fileID,'%12s\n','dump   		dumpstress all mesh/gran/VTK 1000 post/mesh.vtk stress wear drum');

%%%%%%%insert particles%%%%%%%%%%%%%
fprintf(fileID,'%12s\n','run			100000 upto');
fprintf(fileID,'%12s\n','fix		movecad all move/mesh mesh drum rotate origin 0. 0.05 0.05 axis 1. 0. 0. period 1.');

fprintf(fileID,'%12s\n','undump 		dumpstress');
fprintf(fileID,'%12s\n','undump 		dmp');
fprintf(fileID,'%12s\n','dump   		dumpstress all mesh/gran/VTK 100 post/mesh.vtk stress wear drum');
fprintf(fileID,'%12s\n','dump		dmp all custom/vtk 50 post/conveyor.vtk id type type x y z ix iy iz radius c_contactsearch');
fprintf(fileID,'%12s\n','run			200010 upto');

fclose(fileID);
type in.txt






