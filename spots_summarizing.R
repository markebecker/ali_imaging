# 230728
# Making a new big normalized spots file.

library(dplyr)
library(data.table)
library(tidyr)
  
summarize_spots <- function(path, key) {
  filenames <- list.files(path = path, full.names = TRUE)
    
  process_csv <- function(file) {
    data <- read.csv(file, stringsAsFactors = FALSE)
    trimmed <- data[-(1:3),]
    summarized_data <- trimmed %>%
      group_by(FRAME) %>%
      summarize(file = basename(file),
                freq = n(), 
                med.max.gfp = median(as.numeric(MAX_INTENSITY_CH3)), 
                med.max.spytub = median(as.numeric(MAX_INTENSITY_CH1)), 
                med.max.nv = median(as.numeric(MAX_INTENSITY_CH2)))
    return(summarized_data)
  }
  
  # Process all CSVs and combine the results
  result <- lapply(filenames, process_csv)
  combined_result <- do.call(rbind, result)
  combined_result$FRAME <- as.numeric(combined_result$FRAME)
  combined_result$FRAME <- combined_result$FRAME + 1
  
  # Complete frames with missing data
  all_frames <- key %>%
    group_by(file, video.end.frame) %>%
    summarize(FRAME = as.integer(seq_len(first(video.end.frame)))) %>%
    ungroup()
  
  # Prune all_frames to include only files in combined_result
  all_frames <- semi_join(all_frames, combined_result, by = "file")
  
  #Add rows from all_frames that are not present in combined_result
  merged <- left_join(all_frames, combined_result, by = c('file', 'FRAME'))
  
  # Fill rows with information from key
  filled <- left_join(merged, key, by = "file")
  fillt <- filled %>%
    mutate(freq = replace_na(freq, 0))
  return(fillt)
}


spots_key <- read.csv("D:/spots_key.csv")
path <- "D:/spots_good/"
spotsall <- summarize_spots(path, spots_key)

spotsall$norfreq 
norspots <- (spotsall %>% group_by(file) %>% mutate(norfreq = as.numeric(freq)-min(as.numeric(freq))))

# fix the d5d6 timing issue...
forspots <- norspots %>%
  mutate(FRAME = ifelse(expt.name == 'd5d6' & FRAME > 29, FRAME + 8, FRAME))

forspotss <- subset(forspots, (expt.name != 'd5d6') | (FRAME < 59))
forspotss$time <- (forspotss$FRAME-1)*2+1

fwrite(forspotss, 'D:/230729_goodspots.csv')

norspots <- forspotss
good.spots.summary <- forspotss %>%
  group_by(file) %>%
  summarise(good.peak.spots = max(norfreq), good.peak.spots.time = time[which.max(norfreq)])
good.spots.summary.j <- inner_join(good.spots.summary, spots_key)

fwrite(good.spots.summary.j, 'D:/230729_goodspots_summary.csv')
