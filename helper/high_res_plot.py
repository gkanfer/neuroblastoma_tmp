import matplotlib.pyplot as plt
import spatialdata as sd
import spatialdata_plot

class PunchPlotter:
    def __init__(self, sdata,element_name,legend_loc = 'right margin', dpi=100, axis_off=True):
        self.sdata = sdata
        self.element_name = element_name
        self.sdata_crop = self.crop_sdata()
        self.legend_loc = legend_loc # 'right margin'/"bottom margin"
        self.dpi = dpi
        self.axis_off = axis_off
        
    def crop_sdata(self):
        extent = sd.get_extent(
            self.sdata[self.element_name],
            coordinate_system="global",
            exact=True,
        )

        xmin, xmax = extent["x"]
        ymin, ymax = extent["y"]

        sdata_cropped = self.sdata.query.bounding_box(
            min_coordinate=[xmin, ymin],
            max_coordinate=[xmax, ymax],
            axes=("x", "y"),
            target_coordinate_system="global",
        )
        return sdata_cropped
    
    def plot_punch(self, colname, colname_map, palette=None, dpi=100):
        self.sdata_crop["table"].obs[colname] = self.sdata_crop["table"].obs["cell_id"].map(colname_map)
        self.sdata_crop["table"].obs[colname] = self.sdata_crop["table"].obs[colname].astype("category")

        plt.style.use("dark_background")
        ax = self.sdata_crop.pl.render_images(element="morphology_focus", channel='DAPI',cmap= 'grey', colorbar=False).pl.render_labels("cell_labels", color=colname, palette=palette).pl.show(coordinate_systems="global",title=f" ",figsize=(10, 10),dpi = self.dpi,na_in_legend=False,legend_loc = self.legend_loc,return_ax=True)
        if self.axis_off:
            ax.set_axis_off()
        plt.show()
        
