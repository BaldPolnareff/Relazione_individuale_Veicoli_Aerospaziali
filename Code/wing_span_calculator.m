function [Wing_span] = wing_span_calculator (Aspect_Ratio, Wing_Surface)

    Wing_span = sqrt(Aspect_Ratio * Wing_Surface);

end