function D = GeneralizedElasticMatrix(E,mu,G,t,h,thick,type)

switch thick
    case 'plate'    
        switch type
            case 'Bend'
                D = [   -(E*t)/(mu^2 - 1),                           0, -(E*mu*t)/(mu^2 - 1),                           0,                  0,                     0;
                                        0,    -(E*t^3)/(12*(mu^2 - 1)),                    0, -(E*mu*t^3)/(12*(mu^2 - 1)),                  0,                     0;
                     -(E*mu*t)/(mu^2 - 1),                           0,    -(E*t)/(mu^2 - 1),                           0,                  0,                     0;
                                        0, -(E*mu*t^3)/(12*(mu^2 - 1)),                    0,    -(E*t^3)/(12*(mu^2 - 1)),                  0,                     0;
                                        0,                           0,                    0,                           0, (E*t)/(2*(mu + 1)),                     0;
                                        0,                           0,                    0,                           0,                  0, (E*t^3)/(24*(mu + 1))];
            case 'Shear'
                D = [(5*G*t)/6,0;
                    0, (5*G*t)/6];
        end
    case 'stub'
       switch type
            case 'Bend'
                D = [        -(E*(h + t))/(mu^2 - 1),                                     (E*h*(h + t))/(2*(mu^2 - 1)),      -(E*mu*(h + t))/(mu^2 - 1),                                  (E*h*mu*(h + t))/(2*(mu^2 - 1)),                           0,                                                   0;
                        (E*h*(h + t))/(2*(mu^2 - 1)),       - (E*t^3)/(24*(mu^2 - 1)) - (E*(h + t/2)^3)/(3*(mu^2 - 1)), (E*h*mu*(h + t))/(2*(mu^2 - 1)), - (E*mu*t^3)/(24*(mu^2 - 1)) - (E*mu*(h + t/2)^3)/(3*(mu^2 - 1)),                           0,                                                   0;
                          -(E*mu*(h + t))/(mu^2 - 1),                                  (E*h*mu*(h + t))/(2*(mu^2 - 1)),         -(E*(h + t))/(mu^2 - 1),                                     (E*h*(h + t))/(2*(mu^2 - 1)),                           0,                                                   0;
                     (E*h*mu*(h + t))/(2*(mu^2 - 1)), - (E*mu*t^3)/(24*(mu^2 - 1)) - (E*mu*(h + t/2)^3)/(3*(mu^2 - 1)),    (E*h*(h + t))/(2*(mu^2 - 1)),       - (E*t^3)/(24*(mu^2 - 1)) - (E*(h + t/2)^3)/(3*(mu^2 - 1)),                           0,                                                   0;
                                                   0,                                                                0,                               0,                                                                0,    (E*(h + t))/(2*(mu + 1)),                         -(E*h*(h + t))/(4*(mu + 1));
                                                   0,                                                                0,                               0,                                                                0, -(E*h*(h + t))/(4*(mu + 1)), (E*(4*h^3 + 6*h^2*t + 3*h*t^2 + t^3))/(24*(mu + 1))];
            case 'Shear'
%                 
                D = [(5*G*(h + t))/6,0;
                    0, (5*G*(h + t))/6];
       end
       end
end