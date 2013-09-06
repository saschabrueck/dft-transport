disp('Checking sparse matrix result');

test=readcsv([structure, '-GF.csv']);

for ic=0:NBlock-1

  % GR_{i,i-1}
  if ic ~= 0
    pref=['GR_', num2str(ic), num2str(ic-1)];
    GRiim1_test=full(test(Bmin(ic+1):Bmax(ic+1),Bmin(ic):Bmax(ic)));
    GRiim1_ref=full(squeeze(GRNNp1(ic,:,:)).');
    nzdiff=nonzero_diff(GRiim1_test,GRiim1_ref);
    abs_err=max(max(nzdiff));
    rel_err=max(max(nzdiff ./ GRiim1_ref));
    if abs(abs_err) < precision
      disp(['   ', pref, ': pass     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    else
      disp(['   ', pref, ': fail     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    end
  end

  % GR_{i}
  pref=['GR_', num2str(ic), num2str(ic)];
  GRii_test=full(test(Bmin(ic+1):Bmax(ic+1),Bmin(ic+1):Bmax(ic+1)));
  GRii_ref=full(squeeze(GR(ic+1,:,:)));
  nzdiff=nonzero_diff(GRii_test,GRii_ref);
  abs_err=max(max(nzdiff));
  rel_err=max(max(nzdiff ./ GRii_ref));
  if abs(abs_err) < precision
    disp(['   ', pref, ': pass     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   ', pref, ': fail     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  end

  % GR_{i,i+1}
  if ic ~= NBlock-1
    pref=['GR_', num2str(ic), num2str(ic+1)];
    GRiip1_test=full(test(Bmin(ic+1):Bmax(ic+1),Bmin(ic+2):Bmax(ic+2)));
    GRiip1_ref=full(squeeze(GRNNp1(ic+1,:,:)));
    nzdiff=nonzero_diff(GRiip1_test,GRiip1_ref);
    abs_err=max(max(nzdiff));
    rel_err=max(max(nzdiff ./ GRiip1_ref));
    if abs(abs_err) < precision
      disp(['   ', pref, ': pass     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    else
      disp(['   ', pref, ': fail     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    end
  end
end
