% Note: this script will not run unless file_in and file_out are defined
% Example:
% file_in = '/data/gusev/USERS/rlc47/PROFILE/LOHGIC/LOHGIC/inputs/LOHGIC.input.1.tsv';
% file_out = '/data/gusev/USERS/rlc47/PROFILE/LOHGIC/LOHGIC/outputs/LOHGIC.output.1.tsv';

ddf = 0.005;
ddp = 0.01;

data = dlmread (file_in, '\t', 0, 0);
fid = fopen (file_out, 'wt');
fprintf (fid, 'ploidy\tvaf\tdepth\tpurity\tlohgic somatic weight\tlohgic germline weight\tlohgic model 1\tlohgic weight 1\tlohgic model 2\tlohgic weight 2\n');

for k=1:size(data, 1)

    ploidy = data (k, 1); %Copy-number at variant locus
    freq = data (k, 2); %VAF
    d = data (k, 3); %Total depth
    p0 = data (k, 4); %Purity of sample
    df_ci = data (k, 5); %VAF confidence interval
    dp_ci = data (k, 6); %purity confidence interval

    types = {};
    types{1} = 'Somatic, LOH CN_{mut}=1';
    somatic_idx(1) = 1;
    if (ploidy > 1)
        l = length(types);
        for i=1:ploidy
            types{l+i} = sprintf ('Somatic, CN_{mut}=%i', i);
            if (i > 1); types{l+i} = sprintf ('Somatic LOH, CN_{mut}=%i', i); end
            somatic_idx(l+i) = l+i;
        end
    end
    l = length(types);
    types{l+1} = 'Germline, LOH CN_{mut}=1';
    germline_idx(1) = l+1; 
    if (ploidy > 1)
        l = length(types);
        for i=1:ploidy
            types{l+i} = sprintf ('Germline, CN_{mut}=%i', i);
            if (i > 1); types{l+i} = sprintf ('Germline LOH, CN_{mut}=%i', i); end
            germline_idx(1+i) = l+i; 
        end
    end
    
    [~, df] = binofit(round(d*freq), d, df_ci);
    
    fs = df(1) : ddf : df(2); 
    fs = sort([fs(1:end-1), df(2)]);
    dp = sort(p0*(1-dp_ci) : ddp : p0*(1+dp_ci));
    
    aics = zeros (length(fs), length(dp), size(types, 2));
    
    for j=1:length(dp)
        p = dp(j);
        for k=1:length(fs)
            f = fs(k);
            
            aic = [];
            aic(1) = 2 - 2 * log (binopdf (round(d*f), d, (p)/(2*(1-p)+1*p))); %somatic LOH
            if (ploidy > 1)
                l = length(aic);
                for i=1:ploidy
                    aic(l+i) = 2 - 2 * log (binopdf (round(d*f), d, (i*p)/(2*(1-p)+ploidy*p))); %somatic LOH high CN
                end
            end
            l = length(aic);
            aic(l+1) = 2 - 2 * log (binopdf (round(d*f), d, (1-p+p)/(2*(1-p)+1*p))); %germline LOH high CN;
            if (ploidy > 1)
                l = length(aic);
                for i=1:ploidy
                    aic(l+i) = 2 - 2 * log (binopdf (round(d*f), d, (1-p+i*p)/(2*(1-p)+ploidy*p))); %germline LOH high CN;
                end
            end
    
            w = zeros(1, size(types, 2));
            if (j == 1); ws = zeros(length(fs), length(dp), size(types, 2)); end;
            not_inf = ~isinf(aic);
    
            D = sum(exp(-0.5*(aic(not_inf)-min(aic(not_inf)))));
            w(not_inf) = exp(-0.5*(aic(not_inf)-min(aic(not_inf)))) / D;
    
            ws (k, j, 1:size(types, 2)) = w;
            aics (k, j, 1:size(types, 2)) = aic;
            clear aic w w_sorted_i;
        end
    end
    
    not_nan = ~isnan(aics);
    D = sum(exp(-0.5*(aics(not_nan)-min(aics(not_nan)))));
    ww = exp(-0.5*(aics-min(aics(not_nan)))) / D;
    
    sum_ww = zeros (1, size(types, 2));
    for i=1:size(types, 2)
        w = ww(:,:,i); 
        sum_ww (i) = sum(w(~isnan(w)));
        pred_outs{i} = sprintf ('%s\t%2.2e', types{i}, sum_ww (i));
        leg{i} = sprintf ('%s, w = %2.2f', types{i}, sum_ww (i));
    end
    
    aics_p = zeros (length(fs), size(types, 2));
    wf = zeros (length(fs), length(dp), size(types, 2));
    for j=1:length(dp)
        for i=1:size(types, 2)
            aics_p(1:length(fs), i) = aics(:,j,i); 
        end
        
        for k=1:length(fs)
            aic_f = aics_p(k,:);
            not_nan = ~isnan(aic_f);
            if any(not_nan)
                D = sum(exp(-0.5*(aic_f(not_nan)-min(aic_f(not_nan)))));
                wf(k, j, 1:size(types, 2)) = exp(-0.5*(aic_f-min(aic_f(not_nan)))) / D;
            end
        end
    end
        
    [~, i] = sort (sum_ww, 'descend');
    fprintf (fid, '%i\t%2.2f\t%i\t%2.2f\t%2.2e\t%2.2e\t', ploidy, freq, d, p0, sum(sum_ww(somatic_idx)), sum(sum_ww(germline_idx)));
    fprintf (fid, '%s\t%2.2e\t%s\t%2.2e\n', types{i(1)}, sum_ww(i(1)), types{i(2)}, sum_ww(i(2)));

    clear aic w somatic_idx germline_idx;
end
fclose (fid);
