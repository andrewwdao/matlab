classdef DoAVisualisation < handle
    properties
        tx_pos
        rx_pos
        area_size
        bartlettspect
        yconv_dB
        mvdrspatialspect
        ymvdr_dB
        musicspatialspect
        ymusic_dB
        est_aoa_mvdr
        max_idx_mvdr
        max_pow_mvdr
        est_aoa_music
        max_idx_music
        max_pow_music
    end
    
    methods
        function obj = DoAVisualisation(tx_pos, rx_pos, area_size, bartlettspect, yconv_dB, mvdrspatialspect, ymvdr_dB, musicspatialspect, ymusic_dB, est_aoa_mvdr, max_idx_mvdr, max_pow_mvdr, est_aoa_music, max_idx_music, max_pow_music)
            obj.tx_pos = tx_pos;
            obj.rx_pos = rx_pos;
            obj.area_size = area_size;
            obj.bartlettspect = bartlettspect;
            obj.yconv_dB = yconv_dB;
            obj.mvdrspatialspect = mvdrspatialspect;
            obj.ymvdr_dB = ymvdr_dB;
            obj.musicspatialspect = musicspatialspect;
            obj.ymusic_dB = ymusic_dB;
            obj.est_aoa_mvdr = est_aoa_mvdr;
            obj.max_idx_mvdr = max_idx_mvdr;
            obj.max_pow_mvdr = max_pow_mvdr;
            obj.est_aoa_music = est_aoa_music;
            obj.max_idx_music = max_idx_music;
            obj.max_pow_music = max_pow_music;
        end
        function plotMap(obj)
            figure('Name', 'Spatial Spectrum with Map visualisation', 'WindowState', 'maximized'); clf;
            subplot(2,2,1); hold on;
            plot(obj.tx_pos(:,1), obj.tx_pos(:,2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
            text(obj.tx_pos(:,1) + 2, obj.tx_pos(:,2), 'Tx', 'Color', 'red', 'FontSize', 12);
            plot(obj.rx_pos(:,1), obj.rx_pos(:,2), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
            text(obj.rx_pos(:,1) + 2, obj.rx_pos(:,2), 'Rx', 'Color', 'blue', 'FontSize', 12);
            xlim([0 obj.area_size]);
            ylim([0 obj.area_size]);
            xlabel('X Position (m)');
            ylabel('Y Position (m)');
            title('Map with Tx and Rx Positions');
            legend('Tx Position', 'Rx Position');
            grid on;
            hold off;
        end
        function addMarkers(est_angle, pow_array_db, color)
            [est_pow, ~] = maxk(pow_array_db, length(est_angle)); % Get the selected maximum pow_array_db 
            for i = 1:length(est_angle)
                plot(est_angle(i), est_pow(i), 'Marker', 'o', 'MarkerSize', 10, 'Color', color, 'LineWidth', 2);
                text(est_angle(i)-5, est_pow(i)+1, ['DoA:', num2str(est_angle(i)), '°; P:', num2str(est_pow(i)), 'dB'], 'Color', color, 'LineWidth', 2);
            end
        end
        function plotNormalSpectrum(obj)
            subplot(2,2,[3,4]);
            plot(obj.bartlettspect.ScanAngles, obj.yconv_dB, ...
                 obj.mvdrspatialspect.ScanAngles, obj.ymvdr_dB, ...
                 obj.musicspatialspect.ScanAngles, obj.ymusic_dB, ...
                 'LineWidth', 2);
            xlabel('Angle (degrees)');
            ylabel('Power (dB)');
            legend('Conventional', 'MVDR', 'MUSIC', 'AutoUpdate', 'off');
            grid on;
            title('Spatial Spectrum');
            hold on;
            obj.addMarkers(obj.est_aoa, obj.ypow, 'blue');
            for i = 1:length(obj.est_aoa_mvdr)
                plot(obj.mvdrspatialspect.ScanAngles(obj.max_idx_mvdr(i)), obj.ymvdr_dB(obj.max_idx_mvdr(i)), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
                text(obj.mvdrspatialspect.ScanAngles(obj.max_idx_mvdr(i))-5, obj.ymvdr_dB(obj.max_idx_mvdr(i))+1, ...
                    ['DoA:', num2str(obj.est_aoa_mvdr(i)), '°; P:', num2str(obj.max_pow_mvdr(i)), 'dB'], 'Color', 'blue', 'LineWidth', 2);
            end
            for i = 1:length(obj.est_aoa_music)
                plot(obj.musicspatialspect.ScanAngles(obj.max_idx_music(i)), obj.ymusic_dB(obj.max_idx_music(i)), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
                text(obj.musicspatialspect.ScanAngles(obj.max_idx_music(i))-5, obj.ymusic_dB(obj.max_idx_music(i))+1, ...
                    ['DoA:', num2str(obj.est_aoa_music(i)), '°; P:', num2str(obj.max_pow_music(i)), 'dB'], 'Color', 'red', 'LineWidth', 2);
            end
            hold off;
        end
    end
end