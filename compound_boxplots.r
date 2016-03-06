library(ggplot2)  #Import plotting library ggplot2
df <- read.csv("compounds.csv")  # Read dataframe

std_num_frames <- 3
std_P_thr <- 0.01
print(paste("The number of consecute frames is set to", std_num_frames,
"unless otherwise specified."))
print(paste("The P value threshold is set to", std_P_thr,
"unless otherwise specified."))

#Change column names
colnames(df) <- c("row", "dose", "rd_onset_frame", "compound",
                 "concentration", "run_number", "rd_metric",
                 "numfr_consec", "p_thr")

# Loop through the concentrations to create boxplots comparing the efficacy of
# each compound against dose.
for (conc in c(1, 2, 5, 10)){
    plt <- ggplot(df[df$concentration == conc &
        df$numfr_consec == std_num_frames & df$p_thr == std_P_thr,],
        aes(compound, dose))
    plot_title <- sprintf("Concentration = %d mM", conc)
    plt + ggtitle(plot_title) + geom_boxplot(aes(fill = rd_metric)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x="Compound", y="Dose (kGy)") +
        scale_fill_discrete(name="Radiation\nDamage\nMetric")
    plot_filename <- sprintf("../Plots/RP_Comparisons/Conc_%d_dose.pdf", conc)
    ggsave(plot_filename)
}

# Loop through the concentrations to create boxplots comparing the efficacy of
# each compound against frame number.
for (conc in c(1, 2, 5, 10)){
    plt <- ggplot(df[df$concentration == conc &
        df$numfr_consec == std_num_frames & df$p_thr == std_P_thr,],
        aes(compound, rd_onset_frame))
    plot_title <- sprintf("Concentration = %d mM", conc)
    plt + ggtitle(plot_title) + geom_boxplot(aes(fill = rd_metric)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        labs(x="Compound", y="Frame Number") +
        scale_fill_discrete(name="Radiation\nDamage\nMetric")
    plot_filename <- sprintf("../Plots/RP_Comparisons/Conc_%d_frame_num.pdf",
    conc)
    ggsave(plot_filename)
}

# Loop through the compounds to see the effect of concentration on the efficacy
# of that particular compound.
for (cmpd in c("Sucrose", "Trehalose", "TEMPO", "Ascorbate", "Glycerol",
               "Ethylene Glycol", "No Protection", "Sodium Nitrate", "DTT")){
    plt <- ggplot(df[df$compound == cmpd & df$numfr_consec == std_num_frames &
        df$p_thr == std_P_thr,], aes(factor(concentration), dose))
    plot_title <- sprintf("%s", cmpd)
    plt + ggtitle(plot_title) + geom_boxplot(aes(fill = rd_metric)) +
    labs(x="Concentration (mM)", y="Dose (kGy)") +
    scale_fill_discrete(name="Radiation\nDamage\nMetric")
    plot_filename <- sprintf("../Plots/%s/%s_conc_comp.pdf", cmpd, cmpd)
    ggsave(plot_filename)

    # Plot dose against the number of consecutive frames required to be similar
    # before radiation damage is considered to have rendered the dataset
    # too damaged.
    plt <- ggplot(df[df$compound == cmpd &
        df$p_thr == std_P_thr & df$rd_metric == "CorMap",],
        aes(factor(numfr_consec), dose))
    plot_title <- sprintf("%s", cmpd)
    plt + ggtitle(plot_title) +
    geom_boxplot(aes(fill = factor(concentration))) +
    labs(x="Consecutive frames for radiation damage onset threshold",
    y="Dose (kGy)") +
    scale_fill_discrete(name="Concentration\n(mM)")
    plot_filename <- sprintf("../Plots/%s/%s_Num_consec_fr_comp.pdf",
    cmpd, cmpd)
    ggsave(plot_filename)

    # Plot dose against P values Threshold.
    plt <- ggplot(df[df$compound == cmpd &
        df$numfr_consec == std_num_frames & df$rd_metric == "CorMap",],
        aes(factor(p_thr), dose))
    plot_title <- sprintf("%s", cmpd)
    plt + ggtitle(plot_title) +
    geom_boxplot(aes(fill = factor(concentration))) +
    labs(x="P value threshold", y="Dose (kGy)") +
    scale_fill_discrete(name="Concentration\n(mM)")
    plot_filename <- sprintf("../Plots/%s/%s_PThresh_comp.pdf", cmpd, cmpd)
    ggsave(plot_filename)
}
