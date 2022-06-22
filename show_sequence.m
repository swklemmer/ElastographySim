function show_sequence(sonograms, img)
%SHOW_SEQUENCE Shows entire sonogram sequence.
for t = 2:size(sonograms, 1)
    pause(0.05)
    set(img, 'CData', squeeze(sonograms(t, :, :)));
end
end

