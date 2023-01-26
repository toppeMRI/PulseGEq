
t = 0:0.01:5;  % a.u.
c = 1.0;       % a.u.

plot(t, c./(c+t).^2)
xlabel('time (a.u.)')
ylabel('nerve impulse response (a.u.)')
