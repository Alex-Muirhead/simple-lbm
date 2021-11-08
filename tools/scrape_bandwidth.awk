BEGIN { FS = ":" }
/bandwidth/ {
    i += 1
    gsub(/\,/, "")
    bandwidths[i] = $2
}
END {
    printf "\t["
    for(key in bandwidths){
        printf "%.3f, ", bandwidths[key]
    }
    print "]"
}