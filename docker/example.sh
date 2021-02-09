data_folder=`printf 'import apa\nprint(apa.path.data_folder)' | python3`
genomes_folder=`printf 'import pybio\nprint(pybio.path.genomes_folder)' | python3`

cd ~/pybio/genomes
./hg38chr22.download.ensembl98.sh ${genomes_folder}

mkdir ${data_folder}
mkdir ${data_folder}/example
cd ${data_folder}/example
for exp_id in 1 2 3 4 5 6
do
  mkdir e${exp_id}
  wget https://expressrna.org/share/data/example/e${exp_id}/example_e${exp_id}.fastq.gz -O e${exp_id}/example_e${exp_id}.fastq.gz
done

wget https://expressrna.org/share/data/example/example.config -O example.config
wget https://expressrna.org/share/data/example/annotation.tab -O annotation.tab

apa.map.lib -lib_id example
