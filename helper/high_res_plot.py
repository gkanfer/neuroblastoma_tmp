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

        element_names = [self.element_name, "morphology_focus", "cell_labels", "table"]
        sdata_plot = self.sdata.subset(element_names, filter_tables=False)

        sdata_cropped = sdata_plot.query.bounding_box(
            min_coordinate=[xmin, ymin],
            max_coordinate=[xmax, ymax],
            axes=("x", "y"),
            target_coordinate_system="global",
        )
        return sdata_cropped
    
    def plot_punch(
        self,
        colname,
        colname_map,
        palette=None,
        dpi=100,
        active_labels=None,
        inactive_value=float("nan"),
        figsize=(10, 10),
        title=" ",
        legend_loc="default",
        show=True,
        save_path=None,
        save_tight=False,
    ):
        self.sdata_crop["table"].obs[colname] = self.sdata_crop["table"].obs["cell_id"].map(colname_map)
        if active_labels is not None:
            active_labels = list(active_labels)
            is_active = self.sdata_crop["table"].obs[colname].isin(active_labels)
            self.sdata_crop["table"].obs[colname] = self.sdata_crop["table"].obs[colname].where(is_active, inactive_value)
        self.sdata_crop["table"].obs[colname] = self.sdata_crop["table"].obs[colname].astype("category")

        plt.style.use("dark_background")
        if legend_loc == "default":
            legend_loc = self.legend_loc
        groups = None
        render_palette = palette
        if palette is not None:
            groups = active_labels if active_labels is not None else list(palette.keys())
            if isinstance(palette, dict):
                render_palette = [palette[label] for label in groups]
        ax = self.sdata_crop.pl.render_images(element="morphology_focus", channel='DAPI',cmap= 'grey').pl.render_labels("cell_labels", color=colname, groups=groups, palette=render_palette).pl.show(coordinate_systems="global",title=title,figsize=figsize,dpi = self.dpi,na_in_legend=False,legend_loc = legend_loc,colorbar=False,return_ax=True)
        if self.axis_off:
            ax.set_axis_off()
        if save_path is not None:
            if save_tight:
                ax.figure.savefig(save_path, dpi=dpi, bbox_inches="tight", pad_inches=0.02)
            else:
                ax.figure.savefig(save_path, dpi=dpi)
        if show:
            plt.show()
        return ax
        
