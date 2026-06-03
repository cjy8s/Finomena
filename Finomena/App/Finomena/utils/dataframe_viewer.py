import pandas as pd
from PySide6.QtWidgets import QWidget, QVBoxLayout, QTableWidget, QTableWidgetItem
from PySide6.QtCore import Signal

class DataFrameViewerWidget(QWidget):
    """A widget to display a pandas DataFrame in a QTableWidget."""
    data_passed_through = Signal(object)

    def __init__(self, parent=None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self.table = QTableWidget()
        layout.addWidget(self.table)

    def load_data(self, df):
        """Loads a DataFrame into the table and then passes it on."""
        if df is None or df.empty:
            self.table.setRowCount(0)
            self.table.setColumnCount(0)
            print("DataFrameViewerWidget: Received empty dataframe.")
        else:
            # Display only the first 500 rows for performance
            df_to_display = df.head(500)
            self.populate_table(df_to_display)
            print(f"DataFrameViewerWidget: Displaying first {len(df_to_display)} of {len(df)} rows.")

        # Pass the ORIGINAL, complete dataframe to the next step
        self.data_passed_through.emit(df)

    def populate_table(self, df):
        """Fills the QTableWidget with data from the DataFrame."""
        self.table.setRowCount(df.shape[0])
        self.table.setColumnCount(df.shape[1])
        self.table.setHorizontalHeaderLabels(df.columns)

        for row in range(df.shape[0]):
            for col in range(df.shape[1]):
                item = QTableWidgetItem(str(df.iat[row, col]))
                self.table.setItem(row, col, item)
        
        self.table.resizeColumnsToContents()
