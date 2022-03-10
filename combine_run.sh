#Run DEM simulation with uniform PSD
mpirun -np 1 lmp_auto < in.uniform



#Run matlab code to obtain beta_star
matlab -nodisplay -nosplash -nodesktop -r "run('<Working_folder_address>/beta_star_extraction.m');exit;"
mv "0.txt" "DEM_collision_mat_uniform.txt"
mv "1.txt" "DEM_ntotal_uniform.txt"


#Main PBM-DEM coupling framework starts

mpirun -np 1 lmp_auto < in.monoPSD
mv "0.txt" "DEM_collision_mat_loop_0.txt"
mv "1.txt" "DEM_ntotal_loop_0.txt"

flag=1

#loop starts
for i in {1..200}
do
   matlab -nodisplay -nosplash -nodesktop -r "run('<Working_folder_address>/PBM_Main("$i").m');exit;"
   
   mpirun -np 1 lmp_auto < in.txt
   mv "0.txt" "DEM_collision_mat_loop_"$i".txt"
   mv "1.txt" "DEM_ntotal_loop_"$i".txt"
   mv "in.txt" "in.txt"$i""
   mv "log.liggghts" "log.liggghts_"$i""
   cp -r "post" "post_"$i""
   echo "$i loops done"
   val=`cat text.txt`
   if [ $val -eq $flag ] 
   then 
    break
   fi 
done

