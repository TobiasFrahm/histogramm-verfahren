function [out] = rcModell(p, omega)
    Hej = p(1) + p(2)./(1 + 1i.*p(2).*p(3).*omega);
    out = complex(real(Hej), imag(Hej)); 
end