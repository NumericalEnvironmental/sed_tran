# sed_tran
This is a Julia language-based script designed to model sediment transport, delineated by particle size class, through adjacent segments of a one-dimensional river in response to changes in sediment-carrying capacity (calculated using literature empirical expressions for suspended load and bed load transport). Specifically, the model deposits or mobilizes sediments in response to changes in river discharge and water column height caused by storm events and dredging. The river bed is divided into thin horizontal layers within each river section to dynamically track the sediment bed composition and thickness over time.

The following tab-delimited input files are required:

* cells.txt - river segment length and width
* dredge.txt - dredging events (time, depth)
* model_params.txt - miscellaneous model parameters (e.g., time step size, upstream water depth, discharge distribution)
* sediments.txt - sediment properties, by size class
* sources.txt - special, external sediment sources - by size class and by cell; applicable mass fluxes are listed as a function of time
* stratigraphy.txt - initial river bed composition, by cell

More background information can be found here: https://numericalenvironmental.wordpress.com/2016/07/11/sediment-transport/

Email me with questions at walt.mcnab@gmail.com. 

THIS CODE/SOFTWARE IS PROVIDED IN SOURCE OR BINARY FORM "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

