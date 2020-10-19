function a = angle(vec1,vec2)
if any(isnan(vec1)) || any(isnan(vec2))
    a = 2;
    return
end
a = abs(dot(vec1,vec2)/(norm(vec1)*norm(vec2)));
end
