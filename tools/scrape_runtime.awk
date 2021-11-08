BEGIN { FS = ":" }
/clock runtime/ {
    i += 1
    gsub(/\,/, "")
    runtimes[i] = $2
}
END {
    printf "\t["
    for(key in runtimes){
        printf "%.3f, ", runtimes[key]
    }
    print "]"
}