base_dir = Dir.pwd

require "pp"

desc "Run demo"
task :default => [:demo]

desc "Run demo"
task :demo do

  base_cmd  = "R --slave --vanilla -f orenogb.R --args"

  # mm10
  bamdir    = '~/Dropbox_Riken/Public_ACCCBiT/Data/Quartz-Seq/bam'
  bam_files = "#{bamdir}/Quartz_01.bam,#{bamdir}/Quartz_02.bam"

  # hg19
  hg19_bamdir    = '~/Dropbox_Riken/Public_ACCCBiT/Data/Smart-Seq2/bam'
  hg19_bam_files = "#{hg19_bamdir}/Smart-Seq2_01.bam,#{hg19_bamdir}/Smart-Seq2_02.bam"

  cmd = "time #{base_cmd} coordination mm10 chr17 35400000 35600000 1 #{bam_files} demo/demo.pdf"
  sh cmd

  cmd = "time #{base_cmd} coordination mm10 chr17 3.55e7+2880 3.55e7+16079 1 #{bam_files} demo/demo2.pdf"
  sh cmd

  cmd = "time #{base_cmd} coordination mm10 chr17 35502880 35516079 1/200 #{bam_files} demo/demo3.pdf"
  sh cmd

  cmd = "time #{base_cmd} gene mm10 Pou5f1 1 #{bam_files} demo/demo4.pdf"
  sh cmd

  cmd = "time #{base_cmd} gene hg19 CDK2 1 #{hg19_bam_files} demo/demo5.pdf"
  sh cmd    

  Dir.glob("demo/*.pdf") do |file|
  	png = "demo/" + File.basename(file, ".pdf") + ".png"
  	cmd = "sips -s format png #{file} --out #{png}"
  	sh cmd
  end  	
end
