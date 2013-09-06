precision = 1e-5
disp('rGF stage1');
for ic=NBlock-1:-1:0

  im = ic + 1;

  pref = ['gR_', num2str(ic)];
  filename = [pref, '.csv'];
  gR_im_ref=full(squeeze(gR(im,:,:)));
  gR_im_test=full(readcsv(filename));
  diff=gR_im_ref-gR_im_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ gR_im_ref));
  if abs(abs_err) < precision
    disp(['   ', pref, ': pass                                      abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   ', pref, ': fail                                      abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  end

end

disp('rGF stage2');

disp('Iteration 0 (conversion from gR)')
GR_ic_ref=full(squeeze(GR(1,:,:)));
GR_ic_test=readcsv('GR_0.csv');
diff=GR_ic_ref-GR_ic_test;
abs_err=max(max(diff));
rel_err=max(max(diff ./ GR_ic_ref));
if abs(abs_err) < precision
  disp(['   GR_0 (written): pass                                  abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
else
  disp(['   GR_0 (written): fail                                  abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  return;
end

for ic=1:NBlock-1

  disp(['Iteration ', num2str(ic)]);

  T_im1i_ref=full(H(Bmin(ic):Bmax(ic),Bmin(ic+1):Bmax(ic+1)));
  T_iim1_ref=full(H(Bmin(ic+1):Bmax(ic+1),Bmin(ic):Bmax(ic)));
  GR_ic_ref=full(squeeze(GR(ic,:,:)));

  % GR_0
  pref = ['GR_', num2str(ic-1)];
  filename = [pref, '_it', num2str(ic), '.csv'];
  GR_ic_test=full(readcsv(filename)).';
  diff=GR_ic_ref-GR_ic_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ GR_ic_ref));
  if abs(abs_err) < precision
    disp(['   ', pref, ': pass                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   ', pref, ': fail                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    return;
  end

  % gR_1
  pref = ['gR_', num2str(ic)];
  filename = [pref, '_it', num2str(ic), '.csv'];
  gR_im_ref=full(squeeze(gR(ic+1,:,:)));
  gR_im_test=full(readcsv(filename));
  diff=gR_im_ref-gR_im_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ gR_im_ref));
  if abs(abs_err) < precision
    disp(['   ', pref, ': pass                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   ', pref, ': fail                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    return;
  end

  % T_01
  pref = ['T_', num2str(ic-1), num2str(ic)];
  filename = [pref, '.csv'];
  T_im1i_test=zeros(Bmax(ic)-Bmin(ic)+1,Bmax(ic+1)-Bmin(ic+1)+1);
  tmp = full(readcsv(filename));
  tmp_size = size(tmp);
  T_im1i_test(1:tmp_size(1),1:tmp_size(2))=tmp;
  diff=T_im1i_ref-T_im1i_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ T_im1i_ref));
  if abs(abs_err) < precision
    disp(['   ', pref, ': pass                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   ', pref, ': fail                                            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    return;
  end

  % GR_0 * T_01
  A_ref=GR_ic_ref*T_im1i_ref;
  A_test=full(readcsv(['A_', num2str(ic), '.csv']));
  diff=A_ref-A_test.';
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ A_ref));
  if abs(abs_err) < precision
    disp(['   GR_im1 * T_im1i: pass                                 abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   GR_im1 * T_im1i: fail                                 abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % GR_0 * T_01 * gR_1
  B_ref=GR_ic_ref*T_im1i_ref*gR_im_ref;
  B_test=full(readcsv(['B_', num2str(ic), '.csv']));
  diff=B_ref-B_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ B_ref));
  if abs(abs_err) < precision
    disp(['   GR_im1 * T_im1i * gR_i: pass                          abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   GR_im1 * T_im1i * gR_i: fail                          abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % - GR_0 * T_01 * gR_1
  D_ref=-1*GR_ic_ref*T_im1i_ref*gR_im_ref;
  D_test=full(readcsv(['D_', num2str(ic), '.csv']));
  diff=D_ref-D_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ D_ref));
  if abs(abs_err) < precision
    disp(['   GRNNp1 (inside): -GR_0 * T_01 * gR_1: pass            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   GRNNp1 (inside): -GR_0 * T_01 * gR_1: fail            abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % T_10 * GR_0 * T_01 * gR_1
  C_ref=T_iim1_ref*GR_ic_ref*T_im1i_ref*gR_im_ref;
  C_test=full(readcsv(['C_', num2str(ic), '.csv']));
  diff=C_ref-C_test.';
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ C_ref));
  if abs(abs_err) < precision
    disp(['   T_10 * GR_0 * T_01 * gR_1: pass                       abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   T_10 * GR_0 * T_01 * gR_1: fail                       abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % gR_1 * T_10 * GR_0 * T_01 * gR_1
  E_ref=gR_im_ref*T_iim1_ref*GR_ic_ref*T_im1i_ref*gR_im_ref;
  E_test=full(readcsv(['E_', num2str(ic), '.csv']));
  diff=E_ref-E_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ E_ref));
  if abs(abs_err) < precision
    disp(['   gR_1 * T_10 * GR_0 * T_01 * gR_1: pass                abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   gR_1 * T_10 * GR_0 * T_01 * gR_1: fail                abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % gR_1 + gR_1 * T_10 * GR_0 * T_01 * gR_1
  F_ref=gR_im_ref+gR_im_ref*T_iim1_ref*GR_ic_ref*T_im1i_ref*gR_im_ref;
  F_test=full(readcsv(['F_', num2str(ic), '.csv']));
  diff=F_ref-F_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ F_ref));
  if abs(abs_err) < precision
    disp(['   GR: gR_1 + gR_1 * T_10 * GR_0 * T_01 * gR_1: pass     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp(['   GR: gR_1 + gR_1 * T_10 * GR_0 * T_01 * gR_1: fail     abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % GR written
  pref=(['   GR_', num2str(ic)]);
  GR_ref=squeeze(GR(ic+1,:,:));
  GR_wr_test=full(readcsv(['GR_', num2str(ic), '.csv']));
  diff=GR_ref-GR_wr_test;
  abs_err=max(max(diff));
  rel_err=max(max(diff ./ GR_ref));
  if abs(abs_err) < precision
    disp([pref, ' (written): pass                                  abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
  else
    disp([pref, ' (written): fail                                  abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    %return;
  end

  % GRNNp1 written
  if ic ~= NBlock - 1
    pref=(['   GR_', num2str(ic), num2str(ic+1)]);
    GRNNp1_ref=squeeze(GRNNp1(ic+1,:,:));
    GRNNp1_wr_test=full(readcsv(['GRNNp1_', num2str(ic+1), '.csv']));
    diff=GRNNp1_ref-GRNNp1_wr_test;
    abs_err=max(max(diff));
    rel_err=max(max(diff ./ GRNNp1_ref));
    if abs(abs_err) < precision
      disp([pref, ' (written): pass                                 abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
    else
      disp([pref, ' (written): fail                                 abs_err=', num2str(abs_err), '    rel_err=', num2str(abs(rel_err))]);
      %return;
    end
  end

end

sparse_check
