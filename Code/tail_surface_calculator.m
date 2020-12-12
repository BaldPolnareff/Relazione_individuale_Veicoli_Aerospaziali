function [Surface] = tail_surface_calculator (k_c, Volume_coeff, MAC, Wing_Surface, Economy_dim, Type)

    Type = lower(Type);
    Aspect_Ratio = 9;
    Wing_span = sqrt(Aspect_Ratio * Wing_Surface);
    l_optimum = k_c * sqrt(4 * MAC * Wing_Surface * Volume_coeff / (pi * Economy_dim(1)));

    if strcmp(Type, 'horizontal') == 1

        Tail_Surface = Volume_coeff * MAC * Wing_Surface / l_optimum;
        Surface = Tail_Surface;

    elseif strcmp(Type, 'vertical') == 1

        Vert_Surface = Volume_coeff * (Wing_Surface * Wing_span) / l_optimum;
        Surface = Vert_Surface;

    else 

        error('Wrong input: either choose between Horizontal or Vertical');

    end

end