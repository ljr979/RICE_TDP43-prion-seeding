
# =====================
# CONFIGURATION SECTION
# =====================
#important - toggle this to FALSE if you want to run this on a big dataset/a full dataset. This will STOP all of the plotting of every image (this slows down your computer) I would leave it as TRUE when doing test runs to check parameters at the start of each experiment.
plot_all_figures = False 
#-----------------
#start here by looping over all files in one treatment, then we can extend to more treatments later
import os
os.environ["KMP_DUPLICATE_LIB_OK"] = "TRUE"
import numpy as np
import matplotlib.pyplot as plt
import tifffile
from cellpose import models
from cellpose.plot import mask_overlay
import matplotlib
if not plot_all_figures:
    matplotlib.use('Agg') 
from matplotlib import colors
import random
from skimage import exposure

def plot_cellpose_results(img, masks):
    """
    Display raw image, categorical masks, and overlay side by side.
    
    Parameters
    ----------
    img : np.ndarray
        Raw image (2D or 3D with channels).
    masks : np.ndarray
        Cellpose output masks (same xy size as img).
    """

    # pick one channel if image is multichannel
    if img.ndim == 3 and img.shape[-1] > 1:
        disp = img[..., 0]   # choose first channel (adjust as needed)
    else:
        disp = img

    # normalize for display
    disp_norm = exposure.rescale_intensity(disp, out_range=(0,1))

    # random categorical colormap for masks
    cmap = colors.ListedColormap(np.random.rand(256, 3))

    # overlay using cellpose utilities
    overlay = mask_overlay(disp_norm, masks)

    # make 3-panel figure
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    axes[0].imshow(disp_norm, cmap="gray")
    axes[0].set_title("Raw image")
    axes[0].axis("off")

    axes[1].imshow(masks, cmap=cmap, interpolation="nearest")
    axes[1].set_title("Masks (labels)")
    axes[1].axis("off")

    axes[2].imshow(overlay)
    axes[2].set_title("Overlay")
    axes[2].axis("off")

    plt.tight_layout()

    return fig

def run_cellpose(extracted_part, model_used='cyto2'):
    """Actually run cellpose

    Args:
        extracted_part (np.array   ): 8 bt image contrast adjusted etc.
        model_used (str, optional): the cellpose model you want to use. Defaults to 'cyto2'.

    Returns:
        _type_: masks for all cells in that image.
    """
    model = models.Cellpose(gpu=False, model_type=model_used)
    masks, flows, styles, diams = model.eval(
                extracted_part,                         # pass the image array
                diameter=0,                  # let Cellpose estimate
                flow_threshold=0.9,
                cellprob_threshold=-5,
                anisotropy=1.0,
                do_3D=False,                 # dimensionmode=2D
                stitch_threshold=0.0,
                channels=[0,1]              # [nuclei_channel, cyto_channel] -> 
            )
    
    return masks

def adjust_image_for_cellpose(plot_all_figures, normalize_image, contrast_stretch, gamma_adjust, filename, image):
    """#this function adjusts the images by normalising them first, stretching the contrast and brightening the middle brighten cells so that it's easier for cellpose. Then it combines these channels and overlays them for visualisation (composite image). Then turns into an 8-bit image for cellpose and saves this combined, adjusted intensity image and returns it so we can run cellpose.
        

    Args:
        plot_all_figures (_type_): function to plot the overlays
        normalize_image (_type_): function used to normalise intensity across image
        contrast_stretch (_type_): function used to stretch contrst
        gamma_adjust (_type_): function used to stretch gamma
        filename (_type_): name of file to save
        image (_type_): the raw image (array) to stretch 

    Returns:
        _type_: adjusted 8 bit image to pass to cellpose
    """
    
    gfp_norm = normalize_image(image)


                # Optional: apply mild contrast stretch (1–99th percentile)
    gfp_adj = contrast_stretch(gfp_norm, low=1, high=99)


                # Optional: apply gamma adjustment (e.g. 0.8 brightens midtones)
    gfp_adj = gamma_adjust(gfp_adj, gamma=0.8)


                #combine the two channels but account for some dominant pixels by halving the images
    combined = gfp_adj 

    if plot_all_figures:
        rgb = np.zeros((gfp_adj.shape[0], gfp_adj.shape[1], 3)) 
        rgb[..., 1] = gfp_adj
        plt.figure(figsize=(6,6))
        plt.imshow(rgb, origin='lower')
        plt.title("False-colour overlay (red=ch1, green=ch0)- weighted sum of channels")
        plt.axis('off')
        plt.show()   # Green channel
                

    combined_8bit = (combined / combined.max() * 255).astype(np.uint8)

    plt.figure(figsize=(6,6))
    plt.imshow(combined_8bit, cmap='gray', origin='lower')
    plt.title(f'{filename}_adjusted')

    return combined_8bit

def normalize_image(img):
    #normalise the intensity of the images so that we can combine them into one image for cellpose
    img = img.astype(float)
    return (img - img.min()) / (img.max() - img.min() + 1e-8)

def contrast_stretch(img, low=1, high=99):
    """Stretch image contrast by clipping low/high percentiles."""
    p_low, p_high = np.percentile(img, (low, high))
    stretched = np.clip((img - p_low) / (p_high - p_low + 1e-8), 0, 1)
    return stretched

def gamma_adjust(img, gamma=1.0):
    """Apply gamma correction: gamma <1 brightens, >1 darkens."""
    return np.power(img, gamma)

    

if __name__ == '__main__':

    #input folder for all your images from your experiment.
    input_folder = 'other_analysis_resources/cellprofiler_pipelines/aggregate_identification/Example_images/confocal_parental_seeded/s1/'
    results_path = 'results/'
    #what is the common part of the filename ? 
    file_identifier = 'ch00'

    #walks through all folders in your input folder and finds all .tif tiles
    filepaths = [[f'{root}/{name}' for name in files if '.tif' in name] for root, dirs, files in os.walk(f'{input_folder}/')]
    #filters for the specific identifier you used above, and anything that isn't 'mask' (already analysed)
    filepaths = [item for sublist in filepaths for item in sublist if f'{file_identifier}' in item and 'mask' not in item]

    #make this TRUE if you are doing a test run of a small subset of images rather than the whole dataset
    test = False
    #how many images to test in the test run
    num_tests = 1
    #make this false if you don't want to save the cellpose masks overlaid over your grayscale image
    save_figs = False
    #cellpose model to use
    model_used = 'cyto2'
    #if you want to save your masks / test outputs into a different location to the original data, make this true, and set the 'results_path' above to that new location.
    diff_output_folder = False

    if test:
        filepaths = random_sample = random.sample(filepaths, num_tests)



    #loop over all the filepaths
    for image_path in filepaths:
        #replace all annoying slashes
        image_path = image_path.replace('\\', '/')
        image_path = image_path.replace('//', '/')
        print(f'currently analysing {image_path}')
        #find the name of the image you're working on (so we can name the output the same)
        filename = image_path.split('/')[-1].split('.')[0]

        #define the actual path where we will save this mask - i.e. same place as the original image. This line just removes from the path to our image, the name of the image and instead gives the folder where that image lives
        save_path = f'{"/".join(image_path.split("/")[0:(len(image_path.split("/"))-1)])}/'

        if diff_output_folder:
            save_path = f'{results_path}/'
        #read in the image
        image = tifffile.imread(image_path)

        if test:
            # show the image itself (the original, greyscale image)
            fig, axes = plt.subplots()
            plt.imshow(image, cmap="gray")
            plt.title(f'{filename}')

         #now normalise both channels, adjust contrast, brighten mid tones, overlay g + r channels and then save combined image before running cellpose
        combined_8bit = adjust_image_for_cellpose(plot_all_figures, normalize_image, contrast_stretch, gamma_adjust, filename, image)
        # --- run Cellpose (match Fiji settings) ---
        masks = run_cellpose(combined_8bit, model_used)

        #this line plots the masks overlaid over your original image so you can check if it works nicely. let's assume this will only run if youre testing the parameters, otherwise it will be really slow. If you want to save these figures too, make the 'save_figs' == True so it saves them during the test run.
        if test:
            plot_cellpose_results(combined_8bit, masks)
            if save_figs:
                plt.savefig(f'{save_path}{filename}_mask_overlay.png')
            plt.show()
            plt.close

        tifffile.imwrite(f'{save_path}{filename}_mask.tif', masks.astype(np.uint16))
        
        plt.close()