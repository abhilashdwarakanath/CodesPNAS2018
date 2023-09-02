function out = jitter(in,dt);

jitvals = randi(dt,length(in),1);
%jitvals = randi([-(dt/2),dt/2],length(in),1);

for i = 1:length(jitvals)
    r = rand;
    if r > 0.5
        o = in(i)+jitvals(i);
    else
        o = in(i)-jitvals(i);
    end
    if o <= 0 || o > max(in)
        out(i,1) = in(i);
    else
        out(i,1) = o;
    end
end