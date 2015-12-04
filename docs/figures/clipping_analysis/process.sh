pybio.bamclip -bam "~/apa/data.apa/20150514_LUHMGU3/*/*/*.bam" -image clipping_analysis
cp clipping_analysis.png ~/apa/docs/source/figures
rm *.png *.tab *.pdf
cd ~/apa/docs
make html
