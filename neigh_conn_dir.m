% This code helps you to create the connectivity neighbors struct from the standard neighbor struct

labels = {'Fp1' 'AF3' 'F7' 'F3' 'FC1' 'FC5' 'T7' 'C3', ...
'CP1' 'CP5' 'P7' 'P3' 'Pz' 'PO3' 'O1' 'Oz' 'O2' 'PO4', ...
'P4' 'P8' 'CP6' 'CP2' 'C4' 'T8' 'FC6' 'FC2' 'F4' 'F8', ...
'AF4' 'Fp2' 'Fz' 'Cz'};

labels = labels';
neigh_cmb = cell(32,32);

for ch1 = 1:32
    ch1_lbl = labels{ch1};
    neigh_ch1 = neighbours(ch1).neighblabel;
    for ch2 = 1:32
        if ch1 ~= ch2
            neigh_temp = {}; k = 1;
            ch2_lbl = labels{ch2};
            ch_cmb = strcat(ch2_lbl,'_to_',ch1_lbl);
            label_cmb{ch1,ch2} = ch_cmb;
            neigh_ch2 = neighbours(ch2).neighblabel;

            % add neighbors of the sender to ch1
            for n2 = 1 : length(neigh_ch2)
                if ~strcmp(neigh_ch2{n2}, ch1_lbl)
                    ch_cmb = strcat(neigh_ch2{n2},'_to_',ch1_lbl);
                    neigh_temp{k,1} = ch_cmb;
                    k = k+1;
                end
            end
            % add neighbors of the sender to neighbors of ch1
            for n2 = 1 : length(neigh_ch2)
                for n1 = 1 : length(neigh_ch1)
                    if ~strcmp(neigh_ch1{n1}, neigh_ch2{n2})
                        ch_cmb = strcat(neigh_ch2{n2},'_to_',neigh_ch1{n1});
                        neigh_temp{k,1} = ch_cmb;
                        k = k+1;
                    end
                end
            end

            % delete duplicates


            neigh_cmb{ch1, ch2} = sort(neigh_temp);
        end
    end
end


neighb = struct;
k = 1;
for ch2 = 1:32
    for ch1 = 1:32
        if ch1 ~= ch2
            lbl = label_cmb{ch1,ch2};
            neigh = neigh_cmb{ch1,ch2};
            labels2{k,1} = lbl;
            neighb(k).label = lbl;
            neighb(k).neighblabel = neigh;
            k = k+1;
        end
    end
end

            

