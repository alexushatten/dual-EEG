function [chs,d] = neighbours_256chs(ch)

% GOAL: find neighbours of channel ch on the EGI net defined as < 2.5 cm
% distance

make_xyz

    dxyz = bsxfun(@minus,xyz,xyz(ch,:));            % distance from channels in xyz coordinates
    [d,chs]=sort(sqrt(sum(dxyz.^2,2)));             % 3D-distance from channel
    chs = chs(d<25);                                % closest channels < 2.5 cm
%     chs = setdiff(chs,ch);                          % discard channel itself


function make_xyz
% Creat matrix of xyz coordinates
xyz =[64.6113   46.7215  -34.6834
   59.9586   57.2434  -11.5331
   52.3514   65.0798    6.5544
   42.3215   68.5396   21.3162
   30.9423   65.7764   36.4883
   21.0282   56.6094   49.3503
   10.1437   45.9699   58.6773
         0   30.5573   65.5405
  -10.3833   10.3981   70.0731
   54.5054   66.7276  -25.0304
   46.7463   76.4174   -5.4428
   34.8274   80.5436   12.3218
   22.6561   78.9012   27.5364
   11.0656   71.5443   40.9349
         0   59.6760   51.9786
  -10.1437   45.9699   58.6773
  -20.7938   24.5749   63.7218
   36.8406   86.3959  -17.6175
   25.0629   90.4413   -0.3789
   12.4853   89.5812   15.1643
         0   82.8790   29.2558
  -11.0656   71.5443   40.9349
  -21.0282   56.6094   49.3503
  -30.0472   37.4403   54.8845
   12.1043   96.3471  -14.8805
         0   95.1397   -0.0000
  -12.4853   89.5812   15.1643
  -22.6561   78.9012   27.5364
  -30.9423   65.7764   36.4883
  -38.7025   47.5168   43.9131
         0   97.3515  -28.6907
  -12.1043   96.3471  -14.8805
  -25.0629   90.4413   -0.3789
  -34.8274   80.5436   12.3218
  -42.3215   68.5396   21.3162
  -45.9011   55.7609   30.8486
  -36.8406   86.3959  -17.6175
  -46.7463   76.4174   -5.4428
  -52.3514   65.0798    6.5544
  -56.0918   51.9676   17.5075
  -54.5541   40.6257   31.2758
  -48.9443   30.5070   43.8941
  -41.9268   19.6649   54.2361
  -32.0226    6.3770   63.8234
  -16.8828   -9.4637   71.6516
  -54.5054   66.7276  -25.0304
  -59.9586   57.2434  -11.5331
  -64.1068   45.8679    1.0238
  -65.5296   34.1182   15.7144
  -62.6028   23.6307   30.8882
  -56.3811   12.2208   44.1557
  -47.8120    0.8934   55.0595
  -33.4313  -14.8856   66.4401
  -64.6113   46.7215  -34.6834
  -68.3410   37.2276  -17.4556
  -71.6492   25.9855   -1.2739
  -71.4450   14.9077   14.7808
  -67.4635    5.2688   31.0199
  -60.6557   -6.3695   44.6286
  -48.2712  -20.5051   57.2889
  -70.8181   25.6771  -37.6634
  -74.0176   14.5550  -19.3739
  -75.7803    4.4999   -3.1441
  -74.4822   -3.8319   14.2907
  -70.0610  -12.1875   30.8054
  -61.5548  -24.9746   43.6971
  -72.3400   12.4328  -54.0865
  -74.0430    6.0910  -37.0101
  -75.8709   -8.7087  -21.0820
  -76.2910  -15.1533   -3.4544
  -73.8889  -20.7318   14.4751
  -69.3499  -27.8951   28.9162
  -71.9661    4.2267  -73.2525
  -74.2804  -27.7820  -18.9074
  -72.7350  -34.9705    0.5085
  -69.3526  -39.6674   17.8970
  -61.9024  -42.2912   35.1650
  -51.6079  -40.3274   50.3065
  -37.1861  -37.2791   62.5077
  -20.4766  -31.1659   70.4208
         0  -21.3124   74.3864
  -67.6931  -10.5069  -87.5444
  -69.8087  -42.7980  -31.8670
  -69.6791  -45.2138  -13.7672
  -66.8103  -51.3836    5.8547
  -60.0183  -55.3899   23.9277
  -50.1222  -56.0921   40.6851
  -36.8716  -53.1492   55.0688
  -20.6208  -48.5008   64.6042
         0  -38.1084   71.4952
  -62.3390  -25.0171  -97.7141
  -64.2540  -33.4358  -78.8270
  -65.5430  -45.5937  -59.9714
  -66.4376  -49.4356  -44.5736
  -65.5433  -55.3115  -26.9854
  -62.1626  -61.8543   -8.2573
  -55.3974  -67.2940   10.5391
  -46.5238  -69.1344   27.6169
  -34.1386  -67.6898   43.3214
  -18.8910  -62.2752   55.8736
         0  -68.3312   52.4624
  -60.4608  -43.8565  -91.8872
  -59.8577  -51.0869  -75.5490
  -58.4794  -60.6325  -56.4430
  -58.9475  -63.4012  -40.7301
  -56.2864  -70.2257  -22.9360
  -49.0868  -76.9495   -4.4589
  -40.1065  -80.8834   12.8990
  -28.9686  -78.9489   29.2680
  -16.8085  -74.8061   42.8426
  -52.3849  -60.4029  -88.0139
  -50.7190  -66.6682  -72.8518
  -48.2444  -74.6537  -53.4197
  -46.0168  -78.3811  -37.9828
  -40.2068  -85.1088  -19.1651
  -31.2489  -89.4590   -0.8131
  -18.3862  -89.3077   16.0381
  -10.7064  -84.1825   29.0967
         0  -78.7153   40.0488
  -41.6342  -71.3305  -86.7846
  -39.1977  -77.3570  -69.1551
  -34.4550  -85.2781  -53.5705
  -29.4795  -90.3318  -35.5577
  -19.6848  -95.2794  -17.8630
  -10.1001  -95.1397    0.0000
         0  -91.1990   15.8203
   10.7064  -84.1825   29.0967
   16.8085  -74.8061   42.8426
   18.8910  -62.2752   55.8736
   20.6208  -48.5008   64.6042
   20.4766  -31.1659   70.4208
   16.8828   -9.4637   71.6516
  -31.1170  -81.0129  -84.3910
  -24.8179  -85.3853  -69.7358
  -18.0724  -91.8166  -52.9590
  -10.1252  -95.2326  -34.4726
         0  -97.2960  -17.1658
   10.1001  -95.1397    0.0000
   18.3862  -89.3077   16.0381
   28.9686  -78.9489   29.2680
   34.1386  -67.6898   43.3214
   36.8716  -53.1492   55.0688
   37.1861  -37.2791   62.5077
   33.4313  -14.8856   66.4401
  -16.5558  -88.0369  -83.5388
   -8.3077  -88.1122  -69.4859
         0  -93.6988  -52.5043
   10.1252  -95.2326  -34.4726
   19.6848  -95.2794  -17.8630
   31.2489  -89.4590   -0.8131
   40.1065  -80.8834   12.8990
   46.5238  -69.1344   27.6169
   50.1222  -56.0921   40.6851
   51.6079  -40.3274   50.3065
   48.2712  -20.5051   57.2889
    8.3077  -88.1122  -69.4859
   18.0724  -91.8166  -52.9590
   29.4795  -90.3318  -35.5577
   40.2068  -85.1088  -19.1651
   49.0868  -76.9495   -4.4589
   55.3974  -67.2940   10.5391
   60.0183  -55.3899   23.9277
   61.9024  -42.2912   35.1650
   61.5548  -24.9746   43.6971
   16.5558  -88.0369  -83.5388
   24.8179  -85.3853  -69.7358
   34.4550  -85.2781  -53.5705
   46.0168  -78.3811  -37.9828
   56.2864  -70.2257  -22.9360
   62.1626  -61.8543   -8.2573
   66.8103  -51.3836    5.8547
   69.3526  -39.6674   17.8970
   69.3499  -27.8951   28.9162
   31.1170  -81.0129  -84.3910
   39.1977  -77.3570  -69.1551
   48.2444  -74.6537  -53.4197
   58.9475  -63.4012  -40.7301
   65.5433  -55.3115  -26.9854
   69.6791  -45.2138  -13.7672
   72.7350  -34.9705    0.5085
   73.8889  -20.7318   14.4751
   70.0610  -12.1875   30.8054
   60.6557   -6.3695   44.6286
   47.8120    0.8934   55.0595
   32.0226    6.3770   63.8234
   10.3833   10.3981   70.0731
   41.6342  -71.3305  -86.7846
   50.7190  -66.6682  -72.8518
   58.4794  -60.6325  -56.4430
   66.4376  -49.4356  -44.5736
   69.8087  -42.7980  -31.8670
   74.2804  -27.7820  -18.9074
   76.2910  -15.1533   -3.4544
   74.4822   -3.8319   14.2907
   67.4635    5.2688   31.0199
   56.3811   12.2208   44.1557
   41.9268   19.6649   54.2361
   20.7938   24.5749   63.7218
   52.3849  -60.4029  -88.0139
   59.8577  -51.0869  -75.5490
   65.5430  -45.5937  -59.9714
   75.8709   -8.7087  -21.0820
   75.7803    4.4999   -3.1441
   71.4450   14.9077   14.7808
   62.6028   23.6307   30.8882
   48.9443   30.5070   43.8941
   30.0472   37.4403   54.8845
   60.4608  -43.8565  -91.8872
   64.2540  -33.4358  -78.8270
   74.0430    6.0910  -37.0101
   74.0176   14.5550  -19.3739
   71.6492   25.9855   -1.2739
   65.5296   34.1182   15.7144
   54.5541   40.6257   31.2758
   38.7025   47.5168   43.9131
   62.3390  -25.0171  -97.7141
   67.6931  -10.5069  -87.5444
   71.9661    4.2267  -73.2525
   72.3400   12.4328  -54.0865
   70.8181   25.6771  -37.6634
   68.3410   37.2276  -17.4556
   64.1068   45.8679    1.0238
   56.0918   51.9676   17.5075
   45.9011   55.7609   30.8486
   70.2065   25.8362  -56.9276
   65.4118   41.3561  -57.5944
   69.3363   19.8091  -76.1655
   67.5033    7.8181  -92.8406
   61.7785   -7.2991 -104.6254
   61.3522   50.7721  -69.0332
   64.9302   34.9923  -81.9415
   63.7504   23.5394  -98.0881
   60.2478    8.7707 -111.5861
   53.4995   61.7828  -75.0972
   59.4662   47.5057  -88.0883
   58.2091   37.3639 -102.7429
   54.5687   22.5962 -118.2351
   41.0056   69.9277  -80.0967
   48.9688   58.4304  -93.9211
   50.6532   48.1659 -106.1972
  -41.0056   69.9277  -80.0967
  -48.9688   58.4304  -93.9211
  -50.6532   48.1659 -106.1972
  -53.4995   61.7828  -75.0972
  -59.4662   47.5057  -88.0883
  -58.2091   37.3639 -102.7429
  -54.5687   22.5962 -118.2351
  -61.3522   50.7721  -69.0332
  -64.9302   34.9923  -81.9415
  -63.7504   23.5394  -98.0881
  -60.2478    8.7707 -111.5861
  -65.4118   41.3561  -57.5944
  -70.2065   25.8362  -56.9276
  -69.3363   19.8091  -76.1655
  -67.5033    7.8181  -92.8406
  -61.7785   -7.2991 -104.6254
         0   -3.7908   73.2055];
end
end