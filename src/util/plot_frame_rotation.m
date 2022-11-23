% Plot frame rotation. R0 (optional) rotation of base frame relative to
% world coordinate system, R1 the rotation relative to the base frame
function plot_frame_rotation(R1, R0)
    if nargin==1
        R0 = eye(3);
    end
    O = zeros(3,1); % Origin
    x = R0*[1;0;0];
    y = R0*[0;1;0];
    z = R0*[0;0;1];

    x_new = R1*x;
    y_new = R1*y;
    z_new = R1*z;
    
    % Debug
%     x_new'*y_new
%     x_new'*z_new
%     y_new'*z_new

    figure; hold on;
    colors = ['r'; 'g'; 'b'];
    for i=1:3
        %plot3([O(i), x(i)], [O(i), y(i)], [O(i), z(i)], "Color", colors(i), "LineStyle","--");
        quiver3(O(i), O(i), O(i), x(i), y(i), z(i), "Color", colors(i), "LineStyle",":", "LineWidth", 2);
    end
    for i=1:3
        %plot3([O(i), x_new(i)], [O(i), y_new(i)], [O(i), z_new(i)], "Color", colors(i));
        quiver3(O(i), O(i), O(i), x_new(i), y_new(i), z_new(i), "Color", colors(i), "LineStyle","-", "LineWidth", 2);
    end
    title("Frame rotation")
    lg=legend("$X_{base}$", "$Y_{base}$", "$Z_{base}$", "$X^\prime$", "$Y^\prime$", "$Z^\prime$");
    lg.Interpreter = "latex";
    grid on; grid minor;
    axis equal;
end