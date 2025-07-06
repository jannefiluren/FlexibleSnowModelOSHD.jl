using Plots

snowdepth
SWEtmp

p1 = plot(SWEbuffer, title="SWEbuffer");
p2 = plot(snowdepthbuffer, title="snowdepthbuffer");
plot(p1, p2, layout=(2,1), legend=false)
