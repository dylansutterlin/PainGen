/crnl/data : stackage lent
/home/usr/ : juste les script 
mnt/data/ : performant calcul ( pas sauvgarder! PAs mettre V. originales)

### Connection
pas mnt le crnldata/ -plutot sur crnldata/team


### Slurm
run sbatch.sh 

### Ressources CRNL
>Gitlab 
>convert to BIDS
>REMI (NEurostar franco)

### ssh
ssh -X dylan.sutterlinguindon@10.69.168.62

### Paths

run on : /mnt/data/socialhealth/dSutterlin

data on : /crnldata/socialhealth/projects/2024_PainGen/Behavioral

scripts at  : /home/dylan.sutterlinguindon/genPain/scripts

### Data transfer

scp dylan.sutterlinguindon@10.69.168.6/home/dylan.sutterlinguindon

## Formation 2 cluster

cd crnl/projet_commun/mars2024/ pour commandes to run venv on cluster

>create config file in /.ssh folder to make remote symbolic link to a node
	activate labgate command if remote conenction to vpn
>jupyter lab change the port number e.g., 1234 : ok for testing on frontal node

>VScode : remote ssh extension, need config file, to access remote machine


