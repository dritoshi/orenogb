base_dir = Dir.pwd

require "pp"

desc "Run demo"
task :default => [:demo]

desc "Run demo"
task :demo do
  cmd = 'R --slave --vanilla -f orenogb.R --args chr17 35400000 35600000 1 ~/Sources/Quartz_01.th.rmrRNA.bam,~/Sources/Quartz_02.th.rmrRNA.bam demo/demo.pdf'
  sh cmd

  cmd = 'R --slave --vanilla -f orenogb.R --args chr17 3.55e7+2880 3.55e7+16079 1 ~/Sources/Quartz_01.th.rmrRNA.bam,~/Sources/Quartz_02.th.rmrRNA.bam demo/demo2.pdf'
  sh cmd

  cmd = 'R --slave --vanilla -f orenogb.R --args chr17 35502880 35516079 1/200 ~/Sources/Quartz_01.th.rmrRNA.bam,~/Sources/Quartz_02.th.rmrRNA.bam demo/demo3.pdf'
  sh cmd

  Dir.glob("demo/*.pdf") do |file|
  	png = "demo/" + File.basename(file, ".pdf") + ".png"
  	cmd = "sips -s format png #{file} --out #{png}"
  	sh cmd
  end  	
end
