function cf = gaborwavelet(A,phi,sigmadec,f,t)

    cf = A.*exp(-((t-phi).^2./(2*sigmadec^2))).*cos(2*pi*f.*(t-phi));

end