BEGIN { FS = ":" }
/speed/ {
    i += 1
    gsub(/\,/, "")
    speeds[i] = $2
}
END {
    printf "\t["
    for(key in speeds){
        printf "%.3f, ", speeds[key]
    }
    print "]"
}