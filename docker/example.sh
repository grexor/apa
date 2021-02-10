data_folder=`printf 'import apa\nprint(apa.path.data_folder)' | python3`
comps_folder=`printf 'import apa\nprint(apa.path.comps_folder)' | python3`
genomes_folder=`printf 'import pybio\nprint(pybio.path.genomes_folder)' | python3`
polya_folder=`printf 'import apa\nprint(apa.path.polya_folder)' | python3`

cd ~/pybio/genomes
./hg38chr22.download.ensembl98.sh ${genomes_folder}

mkdir ${data_folder}
mkdir ${comps_folder}
mkdir ${polya_folder}
mkdir ${data_folder}/example
mkdir ${comps_folder}/example
cd ${data_folder}/example
for exp_id in 1 2 3 4 5 6
do
  mkdir e${exp_id}
  wget https://expressrna.org/share/data/example/e${exp_id}/example_e${exp_id}.fastq.gz -O e${exp_id}/example_e${exp_id}.fastq.gz
done

wget https://expressrna.org/share/data/example/example.config -O example.config
wget https://expressrna.org/share/data/example/annotation.tab -O annotation.tab
wget https://expressrna.org/share/comps/example/example.config -O ${comps_folder}/example/example.config

apa.map.lib -lib_id example
apa.bed.multi -lib_id example
apa.polya.makeconfig -lib_id example
apa.polya -poly_id example
cp ${polya_folder}/example.bed.gz ${polya_folder}/example_strong.bed.gz
apa.bed.gene_expression -lib_id example
apa.bed.multi -lib_id example -type expression -poly_id example -upstream 10 -downstream 10
apa.comps -comps_id example
