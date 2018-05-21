library(LOTUS)
library(ggplot2)
data(batnets)

ind <- desired_mets <- c('functional complementarity',
                         'web asymmetry',
                         'Fisher alpha', 'mean number of shared partners',
                         'niche overlap',
                         'nestedness',
                         'discrepancy',
                         'ISA')
m <- metcalcs(networks= batnets, indices = ind, network_level = 'higher')

g <- ggplot(m , aes(x = clustering, y = value, color = network)) +
  geom_point()+
  labs(x = 'clustering') +
  geom_smooth(method = lm, se = T)+
  scale_x_continuous(breaks = seq(91, 98, 1))+
  facet_wrap(~ metric, scales = 'free_y')+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
pdf('example_images/scatterplot.pdf')
g
dev.off()

pdf('example_images/lineplot.pdf')
line_plot(input = m, network = 'network', clustering = 'clustering', metric = 'metric', value = 'value', plotname = 'Batnets example')
dev.off()
