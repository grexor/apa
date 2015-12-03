pybio.bamclip -bam "~/apa/data.apa/20150514_LUHMGU3/*/*/*.bam"
cp bamclip.png ~/apa/docs/source/
rm *.png *.tab *.pdf
cd ~/apa/docs
make html
