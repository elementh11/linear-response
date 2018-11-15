
% uniform !!!FRACTIONAL!!! k-mesh generator
clearvars;

nk1 = 3;
nk2 = 1;
nk3 = 1;


kx = linspace(0, 1, nk1);
ky = linspace(0, 1, nk2);
kz = linspace(0, 1, nk3);
klist = zeros(nk1 * nk2 * nk3, 3);
weightlist = ones(nk1 * nk2 * nk3, 1) / (nk1 * nk2 * nk3);

kc = 0;
for k1counter = 1:nk1
    if(nk1>1)
        k1val = [kx(k1counter) 0 0] ;
    else
        k1val = [0 0 0];
    end
    for k2counter = 1:nk2
        if(nk2>1)
            k2val = [0 ky(k2counter) 0] ;
        else
            k2val = [0 0 0];
        end
        for k3counter = 1:nk3%<--- will need to put custom numbers here
            if(nk3>1)
                k3val = [0 0 kz(k3counter)];
            else
                k3val = [0 0 0];
            end
            kc = kc+1;
            kval = k1val + k2val + k3val;
            klist(kc,:) = kval(1:3);
                if((k1counter == 1 || k1counter == nk1) && nk1 ~= 1)
                    weightlist(kc) = weightlist(kc) / 2;
                end
                if((k2counter == 1 || k2counter == nk2) && nk2 ~= 1)
                    weightlist(kc) = weightlist(kc) / 2;
                end
                if((k3counter == 1 || k3counter == nk3) && nk3 ~= 1)
                    weightlist(kc) = weightlist(kc) / 2;
                end
        end
    end
end

datak.weightlist = weightlist;
datak.kpoints = klist;
save('kfile','datak');%<--- will need to put custom numbers here




