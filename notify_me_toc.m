function notify_me_toc(Toc)
% Notify me with time taken
if (Toc/60)<60
    sendolmail('youremail@example.com','Your PC',['Matlab command terminated in ',num2str(round(Toc/60)),' minutes.']);
else
    sendolmail('youremail@example.com','Your PC',['Matlab command terminated in ',num2str(round(Toc/3600)),' hours and ',num2str(round(rem(Toc,3600)/60)),' minutes.']);
end
end