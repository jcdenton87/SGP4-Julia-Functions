function tleformat(longstr1,longstr2)

for (j = 11:16)
    if (longstr1[j] == ' ')
        longstr1 = string(longstr1[1:j-1],"_",longstr1[j+1:end]);
    end
end

if (longstr1[45] != ' ')
    longstr1 = string(longstr1[1:43],longstr1[45],longstr1[45:end]);
end
longstr1 = string(longstr1[1:44],".",longstr1[46:end]);

if (longstr1[8] == ' ')
    longstr1 = string(longstr1[1:7],"U",longstr1[9:end]);
end

if (longstr1[10] == ' ')
    longstr1 = string(longstr1[1:9],".",longstr1[11:end]);
end

for (j = 46:50)
    if (longstr1[j] == ' ')
        longstr1 = string(longstr1[1:j-1],"0",longstr1[j+1:end]);
    end
end
if (longstr1[52] == ' ')
    longstr1 = string(longstr1[1:51],"0",longstr1[53:end]);
end
if (longstr1[54] != ' ')
    longstr1 = string(longstr1[1:52],longstr1[54],longstr1[54:end]);
end

longstr1 = string(longstr1[1:53],".",longstr1[55:end]);

longstr2 = string(longstr2[1:25],".",longstr2[27:end]);

for (j = 27:33)
    if (longstr2[j] == ' ')
        longstr2 = string(longstr2[1:j-1],"0",longstr2[j+1:end]);
    end
end

if (longstr1[63] == ' ')
    longstr1[63] = "0";
end

if ((length(longstr1) < 68) || (longstr1[68] == ' '))
    longstr1 = string(longstr1[1:67],"0");
end

return longstr1,longstr2
end
